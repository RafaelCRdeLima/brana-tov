#include "mcmc.hpp"

#include "constants.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace brana {
namespace {

struct GWRadiusSample {
    double m1 = 0.0;
    double m2 = 0.0;
    double r1 = 0.0;
    double r2 = 0.0;
};

struct NICERRadiusSample {
    double weight = 1.0;
    double mass = 0.0;
    double radius = 0.0;
};

struct MRBranch {
    std::vector<double> mass;
    std::vector<double> radius;
};

MRBranch stable_branch(const Sequence& seq) {
    MRBranch branch;
    if (seq.stars.size() < 3) {
        return branch;
    }

    auto imax_it = std::max_element(seq.stars.begin(), seq.stars.end(), [](const Star& a, const Star& b) {
        return a.M_sun < b.M_sun;
    });
    const auto imax = static_cast<std::size_t>(std::distance(seq.stars.begin(), imax_it));

    std::vector<std::pair<double, double>> mr;
    mr.reserve(imax + 1);
    for (std::size_t i = 0; i <= imax; ++i) {
        mr.emplace_back(seq.stars[i].M_sun, seq.stars[i].R);
    }
    std::sort(mr.begin(), mr.end());

    for (const auto& item : mr) {
        if (!branch.mass.empty() && item.first <= branch.mass.back()) {
            continue;
        }
        branch.mass.push_back(item.first);
        branch.radius.push_back(item.second);
    }
    return branch;
}

std::vector<GWRadiusSample> load_ligo_virgo_samples(const std::string& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("Nao foi possivel abrir posterior LIGO/Virgo: " + path);
    }

    std::vector<GWRadiusSample> samples;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::istringstream is(line);
        double lambda1 = 0.0;
        double lambda2 = 0.0;
        GWRadiusSample sample;
        if (is >> sample.m1 >> sample.m2 >> lambda1 >> lambda2 >> sample.r1 >> sample.r2) {
            if (std::isfinite(sample.m1) && std::isfinite(sample.m2)
                && std::isfinite(sample.r1) && std::isfinite(sample.r2)) {
                samples.push_back(sample);
            }
        }
    }

    if (samples.empty()) {
        throw std::runtime_error("Posterior LIGO/Virgo sem amostras validas: " + path);
    }
    return samples;
}

std::vector<NICERRadiusSample> thin_nicer_samples(const std::vector<NICERRadiusSample>& samples, int max_samples) {
    if (max_samples <= 0 || static_cast<int>(samples.size()) <= max_samples) {
        return samples;
    }

    std::vector<NICERRadiusSample> thinned;
    thinned.reserve(static_cast<std::size_t>(max_samples));
    const double step = static_cast<double>(samples.size() - 1) / static_cast<double>(max_samples - 1);
    for (int i = 0; i < max_samples; ++i) {
        const auto idx = static_cast<std::size_t>(std::round(static_cast<double>(i) * step));
        thinned.push_back(samples[std::min(idx, samples.size() - 1)]);
    }
    return thinned;
}

// radius_first=true  -> columns: R [km], M [Msun], weight   (J0740 format)
// radius_first=false -> columns: weight, M [Msun], R [km]   (J1231 format)
std::vector<NICERRadiusSample> load_nicer_samples(const std::string& path, int max_samples,
                                                   bool radius_first = false) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("Nao foi possivel abrir posterior NICER: " + path);
    }

    std::vector<NICERRadiusSample> samples;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::istringstream is(line);
        std::vector<double> columns;
        double value = 0.0;
        while (is >> value) {
            columns.push_back(value);
        }

        NICERRadiusSample sample;
        if (columns.size() == 2) {
            sample.radius = columns[0];
            sample.mass   = columns[1];
            sample.weight = 1.0;
        } else if (columns.size() >= 3) {
            if (radius_first) {
                // J0740 layout: R [km], M [Msun], weight
                sample.radius = columns[0];
                sample.mass   = columns[1];
                sample.weight = columns[2];
            } else {
                // J1231 layout: weight, M [Msun], R [km]
                sample.weight = columns[0];
                sample.mass   = columns[1];
                sample.radius = columns[2];
            }
        } else {
            continue;
        }

        if (std::isfinite(sample.weight) && std::isfinite(sample.mass) && std::isfinite(sample.radius)
            && sample.weight > 0.0 && sample.mass > 0.0 && sample.radius > 0.0) {
            samples.push_back(sample);
        }
    }

    if (samples.empty()) {
        throw std::runtime_error("Posterior NICER sem amostras validas: " + path);
    }
    return thin_nicer_samples(samples, max_samples);
}

const std::vector<GWRadiusSample>& ligo_virgo_samples(const std::string& path) {
    static const std::vector<GWRadiusSample> samples = load_ligo_virgo_samples(path);
    return samples;
}

const std::vector<NICERRadiusSample>& nicer_j0740_samples(const ObservationConfig& obs) {
    // J0740 file layout: R [km], M [Msun], weight  -> radius_first = true
    static const std::vector<NICERRadiusSample> samples =
        load_nicer_samples(obs.nicer_j0740_samples_path, obs.max_nicer_samples_per_source, true);
    return samples;
}

const std::vector<NICERRadiusSample>& nicer_j1231_samples(const ObservationConfig& obs) {
    // J1231 file layout: weight, M [Msun], R [km]  -> radius_first = false
    static const std::vector<NICERRadiusSample> samples =
        load_nicer_samples(obs.nicer_j1231_samples_path, obs.max_nicer_samples_per_source, false);
    return samples;
}

double logsumexp(const std::vector<double>& values) {
    if (values.empty()) {
        return -std::numeric_limits<double>::infinity();
    }
    const double max_value = *std::max_element(values.begin(), values.end());
    if (!std::isfinite(max_value)) {
        return -std::numeric_limits<double>::infinity();
    }
    double sum = 0.0;
    for (double value : values) {
        sum += std::exp(value - max_value);
    }
    return max_value + std::log(sum);
}

double logaddexp(double a, double b) {
    if (!std::isfinite(a)) {
        return b;
    }
    if (!std::isfinite(b)) {
        return a;
    }
    const double max_value = std::max(a, b);
    return max_value + std::log(std::exp(a - max_value) + std::exp(b - max_value));
}

double radius_at_mass_stable_branch(const MRBranch& branch, double target_mass) {
    if (branch.mass.size() < 2) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (target_mass < branch.mass.front() || target_mass > branch.mass.back()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    auto upper = std::upper_bound(branch.mass.begin(), branch.mass.end(), target_mass);
    if (upper == branch.mass.begin()) {
        return branch.radius.front();
    }
    if (upper == branch.mass.end()) {
        return branch.radius.back();
    }

    const auto i = static_cast<std::size_t>(std::distance(branch.mass.begin(), upper) - 1);
    const double denom = branch.mass[i + 1] - branch.mass[i];
    const double t = denom != 0.0 ? (target_mass - branch.mass[i]) / denom : 0.0;
    return branch.radius[i] + t * (branch.radius[i + 1] - branch.radius[i]);
}

double empirical_ligo_virgo_log_likelihood(const MRBranch& branch, const ObservationConfig& obs) {
    const auto& samples = ligo_virgo_samples(obs.ligo_virgo_samples_path);
    const double h = obs.radius_bandwidth_km;
    if (h <= 0.0) {
        throw std::runtime_error("radius_bandwidth_km deve ser positivo.");
    }

    std::vector<double> log_terms;
    log_terms.reserve(samples.size());
    int valid_count = 0;

    for (const auto& sample : samples) {
        const double r1_model = radius_at_mass_stable_branch(branch, sample.m1);
        const double r2_model = radius_at_mass_stable_branch(branch, sample.m2);
        if (!std::isfinite(r1_model) || !std::isfinite(r2_model)) {
            continue;
        }
        ++valid_count;

        const double z1 = (r1_model - sample.r1) / h;
        const double z2 = (r2_model - sample.r2) / h;
        log_terms.push_back(-0.5 * (z1 * z1 + z2 * z2));
    }

    if (log_terms.empty()) {
        return -std::numeric_limits<double>::infinity();
    }
    // NB: truncating toward zero reproduces the original integer expression
    // `log_terms.size() < samples.size() / 10` bit-for-bit at fraction = 0.1.
    if (obs.use_support_cut
        && log_terms.size() < static_cast<std::size_t>(
               static_cast<double>(samples.size()) * obs.support_cut_fraction)) {
        return -std::numeric_limits<double>::infinity();
    }

    const double coverage = static_cast<double>(valid_count) / static_cast<double>(samples.size());
    return logsumexp(log_terms) - std::log(static_cast<double>(samples.size()))
        + obs.support_coverage_weight * std::log(std::max(coverage, obs.support_coverage_floor));
}

double empirical_nicer_log_likelihood(
    const MRBranch& branch,
    const std::vector<NICERRadiusSample>& samples,
    const ObservationConfig& obs
) {
    const double bandwidth_km = obs.nicer_radius_bandwidth_km;
    const double support_coverage_weight = obs.support_coverage_weight;
    if (bandwidth_km <= 0.0) {
        throw std::runtime_error("nicer_radius_bandwidth_km deve ser positivo.");
    }

    std::vector<double> log_terms;
    log_terms.reserve(samples.size());
    double log_weight_norm = -std::numeric_limits<double>::infinity();
    double log_valid_weight_norm = -std::numeric_limits<double>::infinity();

    for (const auto& sample : samples) {
        log_weight_norm = logaddexp(log_weight_norm, std::log(sample.weight));
        const double r_model = radius_at_mass_stable_branch(branch, sample.mass);
        if (!std::isfinite(r_model)) {
            continue;
        }
        log_valid_weight_norm = logaddexp(log_valid_weight_norm, std::log(sample.weight));

        const double z = (r_model - sample.radius) / bandwidth_km;
        log_terms.push_back(std::log(sample.weight) - 0.5 * z * z);
    }

    if (log_terms.empty()) {
        return -std::numeric_limits<double>::infinity();
    }
    if (obs.use_support_cut
        && log_terms.size() < static_cast<std::size_t>(
               static_cast<double>(samples.size()) * obs.support_cut_fraction)) {
        return -std::numeric_limits<double>::infinity();
    }

    const double coverage = std::exp(log_valid_weight_norm - log_weight_norm);
    return logsumexp(log_terms) - log_weight_norm
        + support_coverage_weight * std::log(std::max(coverage, obs.support_coverage_floor));
}

double constraint_log_penalty(const Diagnostics& diag, const ObservationConfig& obs) {
    double penalty = 0.0;
    if (!obs.use_constraint_penalties) {
        return penalty;
    }
    if (diag.Mmax < obs.Mmax_min) {
        const double z = (diag.Mmax - obs.Mmax_min) / obs.Mmax_sigma_below;
        penalty += -0.5 * z * z;
    }
    if (diag.Cmax > obs.Cmax_max) {
        const double z = (diag.Cmax - obs.Cmax_max) / obs.Cmax_sigma_above;
        penalty += -0.5 * z * z;
    }
    return penalty;
}

std::vector<MCMCSample> run_single_ensemble(
    const EOS& eos,
    const std::vector<double>& rho_c_values,
    const ObservationConfig& obs,
    const MCMCConfig& config,
    int ensemble,
    std::atomic<int>* completed_steps
) {
    constexpr int ndim = 3;
    std::mt19937 rng(config.seed + 7919u * static_cast<unsigned int>(ensemble));
    std::normal_distribution<double> normal(0.0, 1.0);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    std::vector<std::array<double, 3>> walkers(static_cast<std::size_t>(config.n_walkers));
    std::vector<Diagnostics> diagnostics(static_cast<std::size_t>(config.n_walkers));
    std::vector<double> logp(static_cast<std::size_t>(config.n_walkers), -std::numeric_limits<double>::infinity());

    for (int i = 0; i < config.n_walkers; ++i) {
        for (int tries = 0; tries < 10000; ++tries) {
            std::array<double, 3> theta = config.theta0;
            for (int j = 0; j < 3; ++j) {
                theta[j] += config.initial_width[j] * normal(rng);
            }
            if (obs.fix_alpha_U) {
                theta[1] = obs.fixed_alpha_U;
            }
            if (obs.fix_w_U) {
                theta[2] = obs.fixed_w_U;
            }
            Diagnostics diag;
            const double lp = log_posterior(theta, eos, rho_c_values, obs, &diag);
            if (std::isfinite(lp)) {
                walkers[static_cast<std::size_t>(i)] = theta;
                diagnostics[static_cast<std::size_t>(i)] = diag;
                logp[static_cast<std::size_t>(i)] = lp;
                break;
            }
        }
        if (!std::isfinite(logp[static_cast<std::size_t>(i)])) {
            throw std::runtime_error("Nao foi possivel inicializar todos os walkers em pontos validos.");
        }
    }

    std::vector<MCMCSample> samples;
    samples.reserve(static_cast<std::size_t>(config.n_steps) * static_cast<std::size_t>(config.n_walkers));

    for (int i = 0; i < config.n_steps; ++i) {
        for (int k = 0; k < config.n_walkers; ++k) {
            int partner = k;
            while (partner == k) {
                partner = static_cast<int>(uniform(rng) * config.n_walkers);
                if (partner >= config.n_walkers) {
                    partner = config.n_walkers - 1;
                }
            }

            const double u = uniform(rng);
            const double sqrt_z = (config.stretch_a - 1.0) * u + 1.0;
            const double z = sqrt_z * sqrt_z / config.stretch_a;

            std::array<double, 3> proposal{};
            for (int j = 0; j < 3; ++j) {
                proposal[j] = walkers[static_cast<std::size_t>(partner)][j]
                    + z * (walkers[static_cast<std::size_t>(k)][j] - walkers[static_cast<std::size_t>(partner)][j]);
            }
            if (obs.fix_alpha_U) {
                proposal[1] = obs.fixed_alpha_U;
            }
            if (obs.fix_w_U) {
                proposal[2] = obs.fixed_w_U;
            }

            Diagnostics diag_prop;
            const double logp_prop = log_posterior(proposal, eos, rho_c_values, obs, &diag_prop);
            bool accepted = false;
            const double log_accept = static_cast<double>(ndim - 1) * std::log(z)
                + logp_prop - logp[static_cast<std::size_t>(k)];

            if (std::log(uniform(rng)) < log_accept) {
                walkers[static_cast<std::size_t>(k)] = proposal;
                logp[static_cast<std::size_t>(k)] = logp_prop;
                diagnostics[static_cast<std::size_t>(k)] = diag_prop;
                accepted = true;
            }

            const auto& theta = walkers[static_cast<std::size_t>(k)];
            MCMCSample sample;
            sample.ensemble = ensemble;
            sample.step = i;
            sample.walker = k;
            sample.log10_lambda = theta[0];
            sample.lambda = std::pow(10.0, theta[0]);
            sample.alpha_U = theta[1];
            sample.w_U = theta[2];
            sample.logposterior = logp[static_cast<std::size_t>(k)];
            sample.diag = diagnostics[static_cast<std::size_t>(k)];
            sample.accepted = accepted;
            samples.push_back(sample);
        }

        if (completed_steps) {
            ++(*completed_steps);
        }
    }

    return samples;
}

} // namespace

double log_prior(const std::array<double, 3>& theta, const ObservationConfig& obs) {
    const double log10_lambda = theta[0];
    const double alpha_U = theta[1];
    const double w_U = theta[2];
    constexpr double fixed_tolerance = 1.0e-12;

    if (!(obs.log10_lambda_min <= log10_lambda && log10_lambda <= obs.log10_lambda_max)) {
        return -std::numeric_limits<double>::infinity();
    }
    if (obs.fix_alpha_U && std::abs(alpha_U - obs.fixed_alpha_U) > fixed_tolerance) {
        return -std::numeric_limits<double>::infinity();
    }
    if (obs.fix_w_U && std::abs(w_U - obs.fixed_w_U) > fixed_tolerance) {
        return -std::numeric_limits<double>::infinity();
    }
    if (!(obs.alpha_U_min <= alpha_U && alpha_U <= obs.alpha_U_max)) {
        return -std::numeric_limits<double>::infinity();
    }
    if (!(obs.w_U_min <= w_U && w_U <= obs.w_U_max)) {
        return -std::numeric_limits<double>::infinity();
    }

    double lp = 0.0;
    if (obs.use_gaussian_log10_lambda_prior) {
        const double z = (log10_lambda - obs.log10_lambda_mu) / obs.log10_lambda_sigma;
        lp += -0.5 * z * z;
    }
    if (obs.use_gaussian_alpha_U_prior) {
        const double z = (alpha_U - obs.alpha_U_mu) / obs.alpha_U_sigma;
        lp += -0.5 * z * z;
    }
    if (obs.use_gaussian_w_U_prior) {
        const double z = (w_U - obs.w_U_mu) / obs.w_U_sigma;
        lp += -0.5 * z * z;
    }
    return lp;
}

// Builds the integrator parameters for a sampled point. The closure travels with
// the ObservationConfig so that a chain is fully determined by its profile.
BraneParams brane_params_from(const std::array<double, 3>& theta,
                              const ObservationConfig& obs) {
    BraneParams params{std::pow(10.0, theta[0]), theta[1], theta[2]};
    params.closure = obs.closure;
    params.rho0_geom = obs.closure_rho0_cgs > 0.0
        ? rho_mass_cgs_to_geom_km2(obs.closure_rho0_cgs) : 0.0;
    params.n_suppress = obs.closure_n_suppress;
    return params;
}

// Grid resolution for the monotonicity check. It decides the discarded prior
// fraction, which is a reported result, so it has to be convergence-tested rather
// than assumed. See docs/PHASE_C4.md.
int monotonicity_grid_points() {
    static const int n = [] {
        const char* v = std::getenv("BRANA_MONO_GRID");
        return v != nullptr ? std::max(64, std::atoi(v)) : 512;
    }();
    return n;
}

// True if d p_tot/d rho_tot goes negative anywhere the star can sample. Checked on
// the tabulated grid before integrating, so a rejected point costs a few hundred
// evaluations of cs2_eff instead of a full mass-radius sequence.
bool violates_monotonicity(const EOS& eos, const BraneParams& params,
                           const std::vector<double>& rho_c_values) {
    if (rho_c_values.empty()) {
        return false;
    }
    // rho_min() is 0 once the analytic atmosphere is on, so the lower end comes from
    // the sequence grid instead: rho_c_values spans the same tabulated range the star
    // samples, and below its floor the V1 suppression has driven the Weyl term to
    // ~1e-23 of itself, where no violation is possible.
    const double rho_lo = rho_c_values.front();
    const double rho_hi = rho_c_values.back();
    if (!(rho_hi > rho_lo)) {
        return false;
    }
    const int n_grid = monotonicity_grid_points();
    const double dlog = std::log(rho_hi / rho_lo) / static_cast<double>(n_grid - 1);
    for (int i = 0; i < n_grid; ++i) {
        const double rho = rho_lo * std::exp(dlog * static_cast<double>(i));
        const double p = eos.p_of_rho(rho);
        const double cs2 = eos.cs2_of_rho(rho);
        if (!std::isfinite(cs2)) {
            continue;
        }
        const double cs2_eff = cs2_eff_weyl(rho, p, cs2, params);
        if (!std::isfinite(cs2_eff) || cs2_eff <= 0.0) {
            return true;
        }
    }
    return false;
}

double log_likelihood(
    const std::array<double, 3>& theta,
    const EOS& eos,
    const std::vector<double>& rho_c_values,
    const ObservationConfig& obs,
    Diagnostics* diag_out
) {
    const BraneParams params = brane_params_from(theta, obs);
    if (obs.use_monotonicity_cut && violates_monotonicity(eos, params, rho_c_values)) {
        return -std::numeric_limits<double>::infinity();
    }
    const Sequence seq = mass_radius_sequence(rho_c_values, eos, params, "mcmc", 80.0, 0.08);
    const Diagnostics diag = sequence_diagnostics(seq);
    if (diag_out) {
        *diag_out = diag;
    }

    if (!diag.valid) {
        return -std::numeric_limits<double>::infinity();
    }
    const MRBranch branch = stable_branch(seq);
    if (branch.mass.size() < 2) {
        return -std::numeric_limits<double>::infinity();
    }

    double logL = constraint_log_penalty(diag, obs);

    if (obs.use_ligo_virgo && obs.weight_ligo_virgo != 0.0) {
        const double component = empirical_ligo_virgo_log_likelihood(branch, obs);
        if (!std::isfinite(component)) {
            return -std::numeric_limits<double>::infinity();
        }
        logL += obs.weight_ligo_virgo * component;
    }

    if (obs.use_nicer_j0740 && obs.weight_nicer_j0740 != 0.0) {
        const double component = empirical_nicer_log_likelihood(
            branch,
            nicer_j0740_samples(obs),
            obs
        );
        if (!std::isfinite(component)) {
            return -std::numeric_limits<double>::infinity();
        }
        logL += obs.weight_nicer_j0740 * component;
    }

    if (obs.use_nicer_j1231 && obs.weight_nicer_j1231 != 0.0) {
        const double component = empirical_nicer_log_likelihood(
            branch,
            nicer_j1231_samples(obs),
            obs
        );
        if (!std::isfinite(component)) {
            return -std::numeric_limits<double>::infinity();
        }
        logL += obs.weight_nicer_j1231 * component;
    }

    return logL;
}

LikelihoodBreakdown log_likelihood_breakdown(
    const std::array<double, 3>& theta,
    const EOS& eos,
    const std::vector<double>& rho_c_values,
    const ObservationConfig& obs
) {
    LikelihoodBreakdown out;
    out.prior = log_prior(theta, obs);
    if (!std::isfinite(out.prior)) {
        out.total = -std::numeric_limits<double>::infinity();
        return out;
    }

    const BraneParams params = brane_params_from(theta, obs);
    const Sequence seq = mass_radius_sequence(rho_c_values, eos, params, "breakdown", 80.0, 0.08);
    out.diag = sequence_diagnostics(seq);
    if (!out.diag.valid) {
        out.total = -std::numeric_limits<double>::infinity();
        return out;
    }
    const MRBranch branch = stable_branch(seq);
    if (branch.mass.size() < 2) {
        out.total = -std::numeric_limits<double>::infinity();
        return out;
    }

    out.constraint_penalty = constraint_log_penalty(out.diag, obs);
    out.ligo_virgo = obs.use_ligo_virgo ? empirical_ligo_virgo_log_likelihood(branch, obs) : 0.0;
    out.nicer_j0740 = obs.use_nicer_j0740
        ? empirical_nicer_log_likelihood(branch, nicer_j0740_samples(obs), obs)
        : 0.0;
    out.nicer_j1231 = obs.use_nicer_j1231
        ? empirical_nicer_log_likelihood(branch, nicer_j1231_samples(obs), obs)
        : 0.0;

    out.total = out.prior
        + out.constraint_penalty
        + obs.weight_ligo_virgo * out.ligo_virgo
        + obs.weight_nicer_j0740 * out.nicer_j0740
        + obs.weight_nicer_j1231 * out.nicer_j1231;
    return out;
}

double log_posterior(
    const std::array<double, 3>& theta,
    const EOS& eos,
    const std::vector<double>& rho_c_values,
    const ObservationConfig& obs,
    Diagnostics* diag_out
) {
    const double lp = log_prior(theta, obs);
    if (!std::isfinite(lp)) {
        return -std::numeric_limits<double>::infinity();
    }

    const double ll = log_likelihood(theta, eos, rho_c_values, obs, diag_out);
    if (!std::isfinite(ll)) {
        return -std::numeric_limits<double>::infinity();
    }
    return lp + ll;
}

std::vector<MCMCSample> run_mcmc(
    const EOS& eos,
    const std::vector<double>& rho_c_values,
    const ObservationConfig& obs,
    const MCMCConfig& config
) {
    constexpr int ndim = 3;
    if (config.n_walkers < 2 * ndim) {
        throw std::runtime_error("Sampler tipo emcee exige pelo menos 2*ndim walkers.");
    }
    if (config.n_ensembles < 2) {
        throw std::runtime_error("Gelman-Rubin exige pelo menos dois ensembles independentes.");
    }
    if (config.stretch_a <= 1.0) {
        throw std::runtime_error("Parametro stretch_a deve ser maior que 1.");
    }

    std::vector<std::vector<MCMCSample>> by_ensemble(static_cast<std::size_t>(config.n_ensembles));
    std::atomic<int> completed_steps{0};
    std::atomic<bool> workers_done{false};
    const int total_steps = config.n_steps * config.n_ensembles;

    std::thread progress_thread;
    if (config.show_progress) {
        progress_thread = std::thread([&]() {
            using clock = std::chrono::steady_clock;
            const auto start = clock::now();
            int last_percent = -1;

            while (!workers_done.load()) {
                const int done = completed_steps.load();
                const int percent = total_steps > 0 ? static_cast<int>(100.0 * done / total_steps) : 100;
                if (percent != last_percent) {
                    const auto now = clock::now();
                    const double elapsed = std::chrono::duration<double>(now - start).count();
                    std::cerr << "\rMCMC progress: " << std::setw(3) << percent
                              << "% (" << done << "/" << total_steps
                              << " ensemble-steps, " << std::fixed << std::setprecision(1)
                              << elapsed << " s)" << std::flush;
                    last_percent = percent;
                }
                std::this_thread::sleep_for(std::chrono::milliseconds(250));
            }

            const int done = completed_steps.load();
            const auto now = clock::now();
            const double elapsed = std::chrono::duration<double>(now - start).count();
            std::cerr << "\rMCMC progress: 100% (" << done << "/" << total_steps
                      << " ensemble-steps, " << std::fixed << std::setprecision(1)
                      << elapsed << " s)\n" << std::flush;
        });
    }

#ifdef _OPENMP
    omp_set_num_threads(config.n_threads);
#pragma omp parallel for schedule(dynamic)
#endif
    for (int ensemble = 0; ensemble < config.n_ensembles; ++ensemble) {
        by_ensemble[static_cast<std::size_t>(ensemble)] =
            run_single_ensemble(eos, rho_c_values, obs, config, ensemble, &completed_steps);
    }

    workers_done.store(true);
    if (progress_thread.joinable()) {
        progress_thread.join();
    }

    std::vector<MCMCSample> samples;
    samples.reserve(
        static_cast<std::size_t>(config.n_ensembles)
        * static_cast<std::size_t>(config.n_steps)
        * static_cast<std::size_t>(config.n_walkers)
    );
    // Release each ensemble as it is folded in. Copying while `by_ensemble` was
    // still fully alive held two complete chains at once, which is what set the
    // peak RSS of a production run: 137 MB at 409600 samples, against ~79 MB of
    // resident EoS and observation tables. Moving does not help by itself here
    // (MCMCSample is trivially copyable) -- the saving comes from freeing each
    // source vector immediately, so the transient is one chain plus one
    // ensemble instead of two chains.
    for (auto& ensemble_samples : by_ensemble) {
        samples.insert(samples.end(),
                       std::make_move_iterator(ensemble_samples.begin()),
                       std::make_move_iterator(ensemble_samples.end()));
        std::vector<MCMCSample>().swap(ensemble_samples);
    }
    return samples;
}

void write_mcmc_dat(const std::vector<MCMCSample>& samples, const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Nao foi possivel escrever: " + path);
    }
    out << "# ensemble step walker log10_lambda lambda alpha_U w_U logposterior Mmax Rmax R14 Cmax rhoc_Mmax accepted\n";
    out << std::setprecision(16);
    for (const auto& s : samples) {
        out << s.ensemble << ' ' << s.step << ' ' << s.walker << ' ' << s.log10_lambda << ' ' << s.lambda << ' ' << s.alpha_U << ' ' << s.w_U
            << ' ' << s.logposterior << ' ' << s.diag.Mmax << ' ' << s.diag.Rmax << ' ' << s.diag.R14
            << ' ' << s.diag.Cmax << ' ' << s.diag.rhoc_Mmax << ' ' << (s.accepted ? 1 : 0) << '\n';
    }
}

} // namespace brana
