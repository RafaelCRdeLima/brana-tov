#include "constants.hpp"
#include "eos.hpp"
#include "mcmc.hpp"
#include "tov.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#define BRANA_HAVE_RUSAGE 1
#endif

using namespace brana;

#ifndef BRANA_GIT_COMMIT
#define BRANA_GIT_COMMIT "unknown"
#endif
#ifndef BRANA_GIT_DIRTY
#define BRANA_GIT_DIRTY "unknown"
#endif
#ifndef BRANA_BUILD_FLAGS
#define BRANA_BUILD_FLAGS "unknown"
#endif

namespace {

// FNV-1a over the fully resolved configuration, so that two runs can be
// compared without diffing every field by hand.
std::string config_fingerprint(const std::string& text) {
    std::uint64_t hash = 1469598103934665603ULL;
    for (unsigned char c : text) {
        hash ^= static_cast<std::uint64_t>(c);
        hash *= 1099511628211ULL;
    }
    std::ostringstream os;
    os << std::hex << std::setw(16) << std::setfill('0') << hash;
    return os.str();
}

// FNV-1a over the executable itself. git_commit plus git_dirty does not pin a
// build: every chain of a dirty tree records the same commit, so a rebuild
// mid-campaign would be invisible in the provenance. This makes "same binary?"
// a question the artefact answers on its own.
std::string binary_fingerprint() {
#if defined(__linux__)
    std::ifstream exe("/proc/self/exe", std::ios::binary);
    if (!exe) {
        return "unavailable";
    }
    std::uint64_t hash = 1469598103934665603ULL;
    std::array<char, 65536> buf{};
    while (exe.read(buf.data(), static_cast<std::streamsize>(buf.size())) || exe.gcount() > 0) {
        for (std::streamsize i = 0; i < exe.gcount(); ++i) {
            hash ^= static_cast<std::uint64_t>(static_cast<unsigned char>(buf[static_cast<std::size_t>(i)]));
            hash *= 1099511628211ULL;
        }
    }
    std::ostringstream os;
    os << std::hex << std::setw(16) << std::setfill('0') << hash;
    return os.str();
#else
    return "unavailable";
#endif
}

// Every field of ObservationConfig + MCMCConfig that can change a result.
// Written verbatim into the run directory and hashed into the fingerprint.
std::string effective_config_text(
    const std::string& profile,
    const ObservationConfig& obs,
    const MCMCConfig& config,
    const std::string& eos_path,
    const std::vector<double>& rho_c
) {
    std::ostringstream o;
    o << std::setprecision(17);
    o << "profile " << profile << '\n';
    o << "eos_path " << eos_path << '\n';
    o << "ligo_virgo_samples_path " << obs.ligo_virgo_samples_path << '\n';
    o << "nicer_j0740_samples_path " << obs.nicer_j0740_samples_path << '\n';
    o << "nicer_j1231_samples_path " << obs.nicer_j1231_samples_path << '\n';
    o << "radius_bandwidth_km " << obs.radius_bandwidth_km << '\n';
    o << "nicer_radius_bandwidth_km " << obs.nicer_radius_bandwidth_km << '\n';
    o << "max_nicer_samples_per_source " << obs.max_nicer_samples_per_source << '\n';
    o << "use_ligo_virgo " << obs.use_ligo_virgo << '\n';
    o << "use_nicer_j0740 " << obs.use_nicer_j0740 << '\n';
    o << "use_nicer_j1231 " << obs.use_nicer_j1231 << '\n';
    o << "weight_ligo_virgo " << obs.weight_ligo_virgo << '\n';
    o << "weight_nicer_j0740 " << obs.weight_nicer_j0740 << '\n';
    o << "weight_nicer_j1231 " << obs.weight_nicer_j1231 << '\n';
    o << "support_coverage_weight " << obs.support_coverage_weight << '\n';
    o << "use_support_cut " << obs.use_support_cut << '\n';
    o << "support_cut_fraction " << obs.support_cut_fraction << '\n';
    o << "use_constraint_penalties " << obs.use_constraint_penalties << '\n';
    if (obs.use_monotonicity_cut) {
        o << "use_monotonicity_cut 1" << '\n';
    }
    // Always emitted: the integrator upgrade of Phase C7 is a different numerical
    // model from the archived chains, and the fingerprint should say so.
    {
        const IntegratorConfig& ic = integrator_config();
        o << "integ_rtol " << ic.rtol << '\n';
        o << "integ_atol " << ic.atol << '\n';
        o << "integ_sigma_surface " << ic.sigma_surface << '\n';
        o << "integ_p_stop_frac " << ic.p_stop_frac << '\n';
        o << "integ_rho_join_cgs " << ic.rho_join_cgs << '\n';
        o << "integ_surface_limiter " << (ic.use_surface_limiter ? 1 : 0) << '\n';
    }
    if (obs.support_coverage_floor != 1.0e-12) {
        o << "support_coverage_floor " << obs.support_coverage_floor << '\n';
    }
    // The closure changes the likelihood, so it must reach the fingerprint;
    // otherwise a V2 chain and an Original chain would be indistinguishable.
    // Emitted only when it is not the original closure, so that every chain run
    // before Phase C1 keeps the fingerprint it was recorded with -- adding a
    // line unconditionally would flag config drift where none happened.
    if (obs.closure != WeylClosure::Original) {
        o << "weyl_closure " << static_cast<int>(obs.closure) << '\n';
        o << "closure_rho0_cgs " << obs.closure_rho0_cgs << '\n';
        o << "closure_n_suppress " << obs.closure_n_suppress << '\n';
    }
    o << "log10_lambda_bounds " << obs.log10_lambda_min << ' ' << obs.log10_lambda_max << '\n';
    o << "log10_lambda_gaussian " << obs.use_gaussian_log10_lambda_prior << ' '
      << obs.log10_lambda_mu << ' ' << obs.log10_lambda_sigma << '\n';
    o << "alpha_U_bounds " << obs.alpha_U_min << ' ' << obs.alpha_U_max << '\n';
    o << "alpha_U_gaussian " << obs.use_gaussian_alpha_U_prior << ' '
      << obs.alpha_U_mu << ' ' << obs.alpha_U_sigma << '\n';
    o << "alpha_U_fixed " << obs.fix_alpha_U << ' ' << obs.fixed_alpha_U << '\n';
    o << "w_U_bounds " << obs.w_U_min << ' ' << obs.w_U_max << '\n';
    o << "w_U_gaussian " << obs.use_gaussian_w_U_prior << ' '
      << obs.w_U_mu << ' ' << obs.w_U_sigma << '\n';
    o << "w_U_fixed " << obs.fix_w_U << ' ' << obs.fixed_w_U << '\n';
    o << "Mmax_min " << obs.Mmax_min << ' ' << obs.Mmax_sigma_below << '\n';
    o << "Cmax_max " << obs.Cmax_max << ' ' << obs.Cmax_sigma_above << '\n';
    o << "init_override " << obs.override_init << '\n';
    o << "init_center " << config.theta0[0] << ' ' << config.theta0[1] << ' ' << config.theta0[2] << '\n';
    o << "init_width " << config.initial_width[0] << ' ' << config.initial_width[1]
      << ' ' << config.initial_width[2] << '\n';
    o << "n_steps " << config.n_steps << '\n';
    o << "n_walkers " << config.n_walkers << '\n';
    o << "n_ensembles " << config.n_ensembles << '\n';
    o << "stretch_a " << config.stretch_a << '\n';
    o << "seed " << config.seed << '\n';
    o << "n_rho " << rho_c.size() << '\n';
    o << "rho_c_min_geom " << (rho_c.empty() ? 0.0 : rho_c.front()) << '\n';
    o << "rho_c_max_geom " << (rho_c.empty() ? 0.0 : rho_c.back()) << '\n';
    return o.str();
}

// Written to every chain output directory. Deliberately extends the existing
// diagnostics mechanism instead of adding a parallel one.
void write_run_provenance(
    const std::string& profile,
    const ObservationConfig& obs,
    const MCMCConfig& config,
    const std::string& eos_path,
    const std::vector<double>& rho_c,
    const std::string& output_dir
) {
    const std::string cfg = effective_config_text(profile, obs, config, eos_path, rho_c);
    const std::string fingerprint = config_fingerprint(cfg);

    const auto dir = std::filesystem::path(output_dir);
    {
        std::ofstream out(dir / "effective_config.dat");
        if (!out) {
            throw std::runtime_error("Nao foi possivel escrever effective_config.dat");
        }
        out << "# Fully resolved configuration for this run.\n";
        out << cfg;
    }
    {
        std::ofstream out(dir / "run_provenance.dat");
        if (!out) {
            throw std::runtime_error("Nao foi possivel escrever run_provenance.dat");
        }
        out << "# Provenance for this chain. See docs/PHASE0.md.\n";
        out << "git_commit " << BRANA_GIT_COMMIT << '\n';
        out << "git_dirty " << BRANA_GIT_DIRTY << '\n';
        out << "binary_fingerprint " << binary_fingerprint() << '\n';
        out << "config_fingerprint " << fingerprint << '\n';
        out << "weyl_closure " << static_cast<int>(obs.closure) << '\n';
        out << "closure_rho0_cgs " << obs.closure_rho0_cgs << '\n';
        out << "closure_n_suppress " << obs.closure_n_suppress << '\n';
        out << "compiler " << __VERSION__ << '\n';
        out << "build_flags " << BRANA_BUILD_FLAGS << '\n';
#ifdef _OPENMP
        out << "openmp 1\n";
#else
        out << "openmp 0\n";
#endif
        out << "n_threads " << config.n_threads << '\n';
        out << "profile " << profile << '\n';
        out << "eos_path " << eos_path << '\n';
    }
    std::cout << "provenance: commit=" << BRANA_GIT_COMMIT
              << " dirty=" << BRANA_GIT_DIRTY
              << " config=" << fingerprint << '\n';
}

// Peak resident set size, appended once the run is over (it is not known when the
// rest of the provenance is written). Cheap to record and it is what sizes a
// campaign: the C1 post-processing was OOM-killed at ~6.9 GB on 2026-08-04 with
// no measurement of its own to point at. Not part of the regression comparison,
// which only looks at mcmc_chain.dat.
void append_peak_rss(const std::string& output_dir, double runtime_seconds,
                     std::size_t n_samples, int n_posterior_curves) {
    std::ofstream out(std::filesystem::path(output_dir) / "run_provenance.dat",
                      std::ios::app);
    if (!out) {
        return;  // a missing measurement must never cost a finished chain
    }
    out << "wall_seconds " << std::fixed << std::setprecision(1) << runtime_seconds << '\n';
    // Peak RSS alone does not transfer to a run with a different sample count.
    // What transfers is the ratio, so record the two counts it divides by: the
    // chain held in memory (grows with the run) and the posterior-band draw,
    // which is capped and therefore does NOT grow with it.
    out << "n_samples " << n_samples << '\n';
    out << "n_posterior_curves " << n_posterior_curves << '\n';
#ifdef BRANA_HAVE_RUSAGE
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        // ru_maxrss is in kilobytes on Linux, bytes on macOS.
#if defined(__APPLE__)
        const long peak_kb = usage.ru_maxrss / 1024;
#else
        const long peak_kb = usage.ru_maxrss;
#endif
        out << "peak_rss_kb " << peak_kb << '\n';
        std::cout << "peak RSS = " << peak_kb / 1024.0 << " MB\n";
    }
#else
    out << "peak_rss_kb unavailable\n";
#endif
}


// ---- Phase C3: likelihood server for an external nested sampler -------------
//
// The sampler lives in Python (dynesty); the likelihood stays here. The prior
// transform is served from here too, deliberately: duplicating prior definitions
// on the Python side is exactly the drift that breaks reproducibility, so the
// ObservationConfig of the profile remains the single source of truth.
//
// Protocol, one request per line, one reply per line:
//   T u1 u2 u3   ->  theta1 theta2 theta3          unit cube -> parameters
//   B t1 t2 t3   ->  prior lv j0740 j1231 pen total Mmax R14   (P2 decomposition)
//   L t1 t2 t3   ->  logL Mmax Rmax R14 Cmax rhoc  log-LIKELIHOOD (prior removed)
//   Q            ->  terminates
// logL is -inf outside the prior support, written as "-inf".

// Inverse of the standard normal CDF by bisection on std::erf. A likelihood
// evaluation costs ~30 ms of CPU, so 80 bisection steps are free and this avoids
// depending on an erfinv that <cmath> does not provide.
double normal_quantile(double q) {
    if (q <= 0.0) return -40.0;
    if (q >= 1.0) return 40.0;
    double lo = -40.0, hi = 40.0;
    for (int i = 0; i < 80; ++i) {
        const double mid = 0.5 * (lo + hi);
        const double cdf = 0.5 * (1.0 + std::erf(mid / std::sqrt(2.0)));
        if (cdf < q) { lo = mid; } else { hi = mid; }
    }
    return 0.5 * (lo + hi);
}

// Truncated-Gaussian or uniform inverse CDF on [lo, hi], matching log_prior().
double transform_one(double u, double lo, double hi, bool gaussian, double mu, double sigma,
                     bool fixed, double fixed_value) {
    if (fixed) {
        return fixed_value;
    }
    if (!gaussian) {
        return lo + u * (hi - lo);
    }
    const double a = 0.5 * (1.0 + std::erf((lo - mu) / (sigma * std::sqrt(2.0))));
    const double b = 0.5 * (1.0 + std::erf((hi - mu) / (sigma * std::sqrt(2.0))));
    return mu + sigma * normal_quantile(a + u * (b - a));
}

void run_serve(const EOS& eos, const std::vector<double>& rho_c, const ObservationConfig& obs) {
    std::string line;
    std::cout << std::setprecision(17);
    while (std::getline(std::cin, line)) {
        std::istringstream is(line);
        char cmd = 0;
        if (!(is >> cmd)) {
            continue;
        }
        if (cmd == 'Q') {
            break;
        }
        double a = 0.0, b = 0.0, c = 0.0;
        if (!(is >> a >> b >> c)) {
            // Never stay silent on a malformed request: the client is blocked on
            // readline() and would hang forever. A NumPy 2 caller formatting with
            // repr() sends "np.float64(3.96)" and deadlocked a diagnostic for two
            // hours before this was found.
            std::cout << "ERR bad-request\n";
            std::cout.flush();
            continue;
        }
        if (cmd == 'T') {
            std::cout
                << transform_one(a, obs.log10_lambda_min, obs.log10_lambda_max,
                                 obs.use_gaussian_log10_lambda_prior,
                                 obs.log10_lambda_mu, obs.log10_lambda_sigma, false, 0.0) << ' '
                << transform_one(b, obs.alpha_U_min, obs.alpha_U_max,
                                 obs.use_gaussian_alpha_U_prior,
                                 obs.alpha_U_mu, obs.alpha_U_sigma,
                                 obs.fix_alpha_U, obs.fixed_alpha_U) << ' '
                << transform_one(c, obs.w_U_min, obs.w_U_max,
                                 obs.use_gaussian_w_U_prior,
                                 obs.w_U_mu, obs.w_U_sigma,
                                 obs.fix_w_U, obs.fixed_w_U) << '\n';
        } else if (cmd == 'B') {
            // Per-dataset likelihood decomposition at an arbitrary theta, needed by
            // Phase P2 (residual versus mass). Same terms as likelihood_breakdown.dat.
            const std::array<double, 3> theta{a, b, c};
            const LikelihoodBreakdown bd =
                log_likelihood_breakdown(theta, eos, rho_c, obs);
            std::cout << bd.prior << ' ' << bd.ligo_virgo << ' ' << bd.nicer_j0740
                      << ' ' << bd.nicer_j1231 << ' ' << bd.constraint_penalty
                      << ' ' << bd.total << ' ' << bd.diag.Mmax << ' ' << bd.diag.R14
                      << '\n';
        } else if (cmd == 'L') {
            const std::array<double, 3> theta{a, b, c};
            Diagnostics diag;
            const double lp = log_prior(theta, obs);
            double logl = -std::numeric_limits<double>::infinity();
            if (std::isfinite(lp)) {
                const double post = log_posterior(theta, eos, rho_c, obs, &diag);
                if (std::isfinite(post)) {
                    logl = post - lp;
                }
            }
            std::cout << logl << ' ' << diag.Mmax << ' ' << diag.Rmax << ' '
                      << diag.R14 << ' ' << diag.Cmax << ' ' << diag.rhoc_Mmax << '\n';
        }
        std::cout.flush();
    }
}

void print_usage(const char* exe) {
    std::cout << "Uso: " << exe << " [sequence|sequence-custom|mcmc|breakdown|all] [args]\n"
              << "  sequence:        gera curvas massa-raio para GR e casos de brana/Weyl\n"
              << "                   args: [eos.txt] [output_dir]\n"
              << "  sequence-custom: gera uma curva modificada com lambda alpha_U w_U [saida.dat] [eos.txt] [n_rho] [dr] [closure] [rho0_cgs] [n]\n"
              << "                   exemplo: " << exe << " sequence-custom 1e3 -0.25 0.333 output/data/custom.dat data/SLy.txt 600 0.02\n"
              << "  mcmc:            roda ensembles independentes tipo emcee e bandas posteriores\n"
              << "                   args: [n_steps] [n_walkers] [n_ensembles] [n_threads] [profile] [eos.txt] [output_dir]\n"
              << "  breakdown:       recalcula likelihood_breakdown.dat a partir da cadeia existente\n"
              << "                   args: [profile] [n_steps] [eos.txt] [output_dir]\n"
              << "  reband:          regera banda posterior M-R a partir da cadeia salva\n"
              << "                   args: [chain.dat] [n_sequencias=1000] [eos.txt] [output_dir]\n"
              << "  profile:         dump do perfil radial de uma estrela\n"
              << "                   args: lambda alpha_U w_U rho_c_cgs saida.dat [eos.txt] [dr] [closure] [rho0_cgs] [n]\n"
              << "  all:             roda sequence e mcmc\n";
}

ObservationConfig observation_profile(const std::string& profile) {
    ObservationConfig obs;

    if (profile == "wide-soft" || profile == "no-j1231" || profile == "weighted" || profile == "organic" || profile == "lambda-only" || profile == "perturbative-gr-base") {
        obs.log10_lambda_min = -1.0;
        obs.log10_lambda_max = 8.0;
        obs.alpha_U_min = -4.0;
        obs.alpha_U_max = 2.0;
        obs.w_U_min = -4.0;
        obs.w_U_max = 2.0;
        obs.radius_bandwidth_km = 1.0;
        obs.nicer_radius_bandwidth_km = 1.0;
    }

    if (profile == "no-j1231") {
        obs.use_nicer_j1231 = false;
        obs.weight_nicer_j1231 = 0.0;
    }

    if (profile == "weighted") {
        obs.weight_ligo_virgo = 1.0;
        obs.weight_nicer_j0740 = 0.5;
        obs.weight_nicer_j1231 = 0.5;
    }

    if (profile == "organic") {
        obs.radius_bandwidth_km = 1.2;
        obs.nicer_radius_bandwidth_km = 0.85;
        obs.weight_ligo_virgo = 0.30;
        obs.weight_nicer_j0740 = 1.50;
        obs.weight_nicer_j1231 = 1.25;
        obs.Mmax_min = 2.12;
        obs.Mmax_sigma_below = 0.06;
        obs.Cmax_max = 0.34;
        obs.Cmax_sigma_above = 0.05;
        obs.support_coverage_weight = 8.0;
        obs.max_nicer_samples_per_source = 12000;
        obs.alpha_U_min = -4.0;
        obs.alpha_U_max = 1.5;
        obs.use_gaussian_alpha_U_prior = true;
        obs.alpha_U_mu = 0.0;
        obs.alpha_U_sigma = 0.20;
        obs.w_U_min = -6.0;
        obs.w_U_max = 3.0;
    }

    if (profile == "lambda-only") {
        obs.radius_bandwidth_km = 1.2;
        obs.nicer_radius_bandwidth_km = 0.85;
        obs.weight_ligo_virgo = 0.30;
        obs.weight_nicer_j0740 = 1.50;
        obs.weight_nicer_j1231 = 1.25;
        obs.Mmax_min = 2.12;
        obs.Mmax_sigma_below = 0.06;
        obs.Cmax_max = 0.34;
        obs.Cmax_sigma_above = 0.05;
        obs.support_coverage_weight = 8.0;
        obs.max_nicer_samples_per_source = 12000;
        obs.fix_alpha_U = true;
        obs.fixed_alpha_U = 0.0;
        obs.alpha_U_min = 0.0;
        obs.alpha_U_max = 0.0;
        obs.use_gaussian_alpha_U_prior = false;
        obs.fix_w_U = true;
        obs.fixed_w_U = 0.0;
        obs.w_U_min = 0.0;
        obs.w_U_max = 0.0;
    }

    if (profile == "perturbative-gr-base") {
        obs.radius_bandwidth_km = 1.2;
        obs.nicer_radius_bandwidth_km = 0.85;
        obs.weight_ligo_virgo = 0.30;
        obs.weight_nicer_j0740 = 1.50;
        obs.weight_nicer_j1231 = 1.25;
        obs.Mmax_min = 2.12;
        obs.Mmax_sigma_below = 0.06;
        obs.Cmax_max = 0.34;
        obs.Cmax_sigma_above = 0.05;
        obs.support_coverage_weight = 8.0;
        obs.max_nicer_samples_per_source = 12000;
        obs.log10_lambda_min = -1.0;
        obs.log10_lambda_max = 9.0;
        obs.use_gaussian_log10_lambda_prior = true;
        obs.log10_lambda_mu = 4.0;
        obs.log10_lambda_sigma = 1.5;
        obs.alpha_U_min = -1.5;
        obs.alpha_U_max = 1.5;
        obs.use_gaussian_alpha_U_prior = true;
        obs.alpha_U_mu = 0.0;
        obs.alpha_U_sigma = 0.25;
        obs.w_U_min = -2.0;
        obs.w_U_max = 2.0;
        obs.use_gaussian_w_U_prior = true;
        obs.w_U_mu = 0.0;
        obs.w_U_sigma = 0.50;
    }

    // ---- Phase 0/1: initialization is part of the profile, and the ablation
    // variants of the production likelihood are named profiles. -------------
    //
    // `perturbative-gr`            : production settings (Phase 1 variant A)
    // `perturbative-gr-narrowinit` : legacy initialization used by the six
    //                                non-SLy chains of 2026-07-17; kept only so
    //                                that those chains remain reproducible.
    // `pgr-abl-{B..F}`             : Phase 1 ablations, see docs/PHASE1.md.
    const bool is_pgr_family = profile == "perturbative-gr"
        || profile == "perturbative-gr-narrowinit"
        || profile.rfind("pgr-abl-", 0) == 0
        || profile == "pgr-v1" || profile == "pgr-v1-nocut" || profile == "pgr-v2";

    if (is_pgr_family) {
        obs = observation_profile("perturbative-gr-base");
        obs.override_init = true;
        obs.init_center = {obs.log10_lambda_mu, obs.alpha_U_mu, obs.w_U_mu};
        // The SLy posterior has separated, but allowed, alpha_U--w_U branches.
        // Start every independent ensemble broad enough to populate both.
        obs.init_width = {0.45, 0.40, 0.80};

        if (profile == "perturbative-gr-narrowinit") {
            obs.init_width = {0.45, 0.10, 0.20};
        }
        // Ablations: κ (support-coverage weight), hard support cut, soft penalties.
        if (profile == "pgr-abl-B") { obs.support_coverage_weight = 0.0; }
        if (profile == "pgr-abl-C") { obs.use_support_cut = false; }
        if (profile == "pgr-abl-D") { obs.support_coverage_weight = 0.0; obs.use_support_cut = false; }
        if (profile == "pgr-abl-E") { obs.use_constraint_penalties = false; }
        if (profile == "pgr-abl-F") {
            obs.support_coverage_weight = 0.0;
            obs.use_support_cut = false;
            obs.use_constraint_penalties = false;
        }
        // Phase C1 closure chains. Production likelihood, only the Weyl closure
        // changes. Needed because the closure enters the TOV solution and hence
        // the likelihood: a V2 posterior cannot be obtained by filtering or
        // reweighting the original-closure chain.
        if (profile == "pgr-v2") {
            obs.closure = WeylClosure::V2PressureTracking;
        }
        // pgr-v1 is production; pgr-v1-nocut exists only to price the monotonicity cut.
        if (profile == "pgr-v1" || profile == "pgr-v1-nocut") {
            obs.closure = WeylClosure::V1Suppressed;
            obs.use_monotonicity_cut = (profile == "pgr-v1");
            obs.closure_rho0_cgs = 2.0 * 2.7e14;   // 2 rho_sat, the C1-preferred scale
            obs.closure_n_suppress = 3.0;
        }
        // Phase 3 kappa scan: pgr-kappa-<value>, e.g. pgr-kappa-2, pgr-kappa-16.
    }

    // Same pattern as pgr-kappa-: pgr-floor-1e-6 keeps the production likelihood
    // and only moves the coverage floor, which is what exposes it as a free knob.
    if (profile.rfind("pgr-floor-", 0) == 0) {
        obs = observation_profile("perturbative-gr");
        obs.support_coverage_floor = std::stod(profile.substr(std::string("pgr-floor-").size()));
    }

    if (profile.rfind("pgr-kappa-", 0) == 0) {
        obs = observation_profile("perturbative-gr");
        obs.support_coverage_weight = std::stod(profile.substr(std::string("pgr-kappa-").size()));
    }

    static const std::vector<std::string> known = {
        "default", "wide-soft", "no-j1231", "weighted", "organic", "lambda-only",
        "perturbative-gr", "perturbative-gr-base", "perturbative-gr-narrowinit",
        "pgr-abl-B", "pgr-abl-C", "pgr-abl-D", "pgr-abl-E", "pgr-abl-F",
        "pgr-v1", "pgr-v1-nocut", "pgr-v2",
    };
    if (std::find(known.begin(), known.end(), profile) == known.end()
        && profile.rfind("pgr-kappa-", 0) != 0
        && profile.rfind("pgr-floor-", 0) != 0) {
        throw std::runtime_error("Perfil desconhecido: " + profile);
    }

    return obs;
}

std::vector<double> central_densities_tabulated(const EOS& eos, int n = 120) {
    const double rho_min = std::max(rho_mass_cgs_to_geom_km2(1.0e14), eos.rho_min() * 1.01);
    const double rho_max = std::min(rho_mass_cgs_to_geom_km2(3.0e16), eos.rho_max() * 0.99);
    return logspace(std::log10(rho_min), std::log10(rho_max), n);
}

// Dump the radial profile of one star, plus the causality/monotonicity
// diagnostics of the effective source. See docs/PHASE_CAUSALITY.md.
void run_profile(
    const EOS& eos,
    double lambda_brane,
    double alpha_U,
    double w_U,
    double rho_c_cgs,
    const std::string& path,
    double dr,
    WeylClosure closure = WeylClosure::Original,
    double rho0_cgs = 0.0,
    double n_suppress = 3.0
) {
    const auto parent = std::filesystem::path(path).parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }

    BraneParams params{lambda_brane, alpha_U, w_U};
    params.closure = closure;
    params.rho0_geom = rho0_cgs > 0.0 ? rho_mass_cgs_to_geom_km2(rho0_cgs) : 0.0;
    params.n_suppress = n_suppress;
    // Phase 2 bracket selector; recorded in the profile header below.
    if (const char* b = std::getenv("BRANA_TIDAL_BETA")) { params.tidal_beta = std::atof(b); }
    const double rho_c = rho_mass_cgs_to_geom_km2(rho_c_cgs);
    std::vector<ProfilePoint> profile;
    const Star star = integrate_star(rho_c, eos, params, 80.0, dr, 1.0e-5, &profile);

    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Nao foi possivel escrever: " + path);
    }
    out << "# radial profile  lambda=" << lambda_brane << " alpha_U=" << alpha_U
        << " w_U=" << w_U << " rho_c_cgs=" << rho_c_cgs
        << " closure=" << static_cast<int>(closure)
        << " rho0_cgs=" << rho0_cgs << " n=" << n_suppress << '\n';
    out << "# R_km " << star.R << " M_sun " << star.M_sun << " M_km " << star.M_km
        << " y_R " << star.y_R << '\n';
    out << "# r_km m_km p_geom rho_geom rho_tot_geom p_r_tot_geom cs2 cs2_eff rho_cgs\n";
    out << std::setprecision(16);
    for (const auto& pt : profile) {
        out << pt.r << ' ' << pt.m << ' ' << pt.p << ' ' << pt.rho << ' '
            << pt.rho_tot << ' ' << pt.p_r_tot << ' ' << pt.cs2 << ' '
            << pt.cs2_eff << ' ' << rho_geom_km2_to_cgs(pt.rho) << '\n';
    }
    std::cout << "perfil: R=" << star.R << " km  M=" << star.M_sun << " Msun  "
              << profile.size() << " pontos -> " << path << '\n';
}

void summarize(const Sequence& seq) {
    const Diagnostics d = sequence_diagnostics(seq);
    std::cout << seq.label << '\n';
    if (!d.valid) {
        std::cout << "  sem modelos validos\n";
        return;
    }
    std::cout << "  Mmax = " << d.Mmax << " M_sun\n"
              << "  Rmax = " << d.Rmax << " km\n"
              << "  R1.4 = " << d.R14 << " km\n"
              << "  Cmax = " << d.Cmax << '\n'
              << "  rhoc(Mmax) = " << d.rhoc_Mmax << " g/cm^3\n";
}

void run_sequences(const EOS& eos, const std::vector<double>& rho_c, const std::string& output_dir) {
    std::filesystem::create_directories(output_dir);

    struct Case {
        std::string filename;
        std::string label;
        BraneParams params;
    };

    const std::vector<Case> cases = {
        {"sequence_gr.dat", "GR", {std::numeric_limits<double>::infinity(), 0.0, 0.0}},
        {"sequence_lambda1e2_alpha-025_w-01.dat", "lambda=1e2 alpha_U=-0.25 w_U=-0.1", {1.0e2, -0.25, -0.1}},
        {"sequence_lambda1e2_alpha-050_w033.dat", "lambda=1e2 alpha_U=-0.50 w_U=1/3", {1.0e2, -0.50, 1.0 / 3.0}},
        {"sequence_lambda1e3_alpha-025_w033.dat", "lambda=1e3 alpha_U=-0.25 w_U=1/3", {1.0e3, -0.25, 1.0 / 3.0}},
    };

    for (const auto& item : cases) {
        std::cout << "Calculando " << item.label << "...\n";
        const Sequence seq = mass_radius_sequence(rho_c, eos, item.params, item.label, 80.0, 0.05);
        write_stable_sequence_dat(seq, (std::filesystem::path(output_dir) / item.filename).string());
        summarize(seq);
    }
}

void run_custom_sequence(
    const EOS& eos,
    const std::vector<double>& rho_c,
    double lambda_brane,
    double alpha_U,
    double w_U,
    const std::string& path,
    double dr,
    WeylClosure closure = WeylClosure::Original,
    double rho0_cgs = 0.0,
    double n_suppress = 3.0
) {
    const auto parent = std::filesystem::path(path).parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }

    BraneParams params{lambda_brane, alpha_U, w_U};
    params.closure = closure;
    params.rho0_geom = rho0_cgs > 0.0 ? rho_mass_cgs_to_geom_km2(rho0_cgs) : 0.0;
    params.n_suppress = n_suppress;
    if (const char* b = std::getenv("BRANA_TIDAL_BETA")) { params.tidal_beta = std::atof(b); }

    std::ostringstream label;
    label << "lambda=" << lambda_brane << " alpha_U=" << alpha_U << " w_U=" << w_U
          << " closure=" << static_cast<int>(closure);

    std::cout << "Calculando " << label.str() << "...\n";
    const Sequence seq = mass_radius_sequence(
        rho_c,
        eos,
        params,
        label.str(),
        80.0,
        dr
    );
    write_stable_sequence_dat(seq, path);
    summarize(seq);
    std::cout << "Sequencia modificada em " << path << '\n';
}

double radius_at_mass(const Sequence& seq, double target_mass) {
    if (seq.stars.size() < 3) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    auto imax_it = std::max_element(seq.stars.begin(), seq.stars.end(), [](const Star& a, const Star& b) {
        return a.M_sun < b.M_sun;
    });
    const auto imax = static_cast<std::size_t>(std::distance(seq.stars.begin(), imax_it));

    std::vector<std::pair<double, double>> mr;
    for (std::size_t i = 0; i <= imax; ++i) {
        mr.emplace_back(seq.stars[i].M_sun, seq.stars[i].R);
    }
    std::sort(mr.begin(), mr.end());
    if (target_mass < mr.front().first || target_mass > mr.back().first) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    for (std::size_t i = 0; i + 1 < mr.size(); ++i) {
        if (target_mass >= mr[i].first && target_mass <= mr[i + 1].first) {
            const double denom = mr[i + 1].first - mr[i].first;
            const double t = denom != 0.0 ? (target_mass - mr[i].first) / denom : 0.0;
            return mr[i].second + t * (mr[i + 1].second - mr[i].second);
        }
    }
    return std::numeric_limits<double>::quiet_NaN();
}

double quantile(std::vector<double> values, double q) {
    values.erase(std::remove_if(values.begin(), values.end(), [](double x) { return !std::isfinite(x); }), values.end());
    if (values.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    std::sort(values.begin(), values.end());
    const double pos = q * static_cast<double>(values.size() - 1);
    const auto lo = static_cast<std::size_t>(std::floor(pos));
    const auto hi = static_cast<std::size_t>(std::ceil(pos));
    const double t = pos - static_cast<double>(lo);
    return values[lo] + t * (values[hi] - values[lo]);
}

double mean(const std::vector<double>& values) {
    if (values.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    double sum = 0.0;
    for (double value : values) {
        sum += value;
    }
    return sum / static_cast<double>(values.size());
}

double sample_value(const MCMCSample& sample, int parameter) {
    if (parameter == 0) {
        return sample.log10_lambda;
    }
    if (parameter == 1) {
        return sample.alpha_U;
    }
    return sample.w_U;
}

double rhat_for_parameter(const std::vector<std::vector<double>>& chains) {
    const std::size_t m = chains.size();
    const std::size_t n = chains.empty() ? 0 : chains.front().size();
    if (m < 2 || n < 2) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    std::vector<double> means;
    std::vector<double> variances;
    means.reserve(m);
    variances.reserve(m);

    for (const auto& chain : chains) {
        const double chain_mean = mean(chain);
        double sumsq = 0.0;
        for (double value : chain) {
            const double delta = value - chain_mean;
            sumsq += delta * delta;
        }
        means.push_back(chain_mean);
        variances.push_back(sumsq / static_cast<double>(n - 1));
    }

    const double mean_of_means = mean(means);
    double between_sum = 0.0;
    for (double chain_mean : means) {
        const double delta = chain_mean - mean_of_means;
        between_sum += delta * delta;
    }

    const double W = mean(variances);
    const double B = static_cast<double>(n) * between_sum / static_cast<double>(m - 1);
    if (W <= 0.0 || !std::isfinite(W)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double var_hat = (static_cast<double>(n - 1) / static_cast<double>(n)) * W + B / static_cast<double>(n);
    return std::sqrt(var_hat / W);
}

double autocorr_time_for_chain(const std::vector<double>& chain) {
    const std::size_t n = chain.size();
    if (n < 3) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double chain_mean = mean(chain);
    double var = 0.0;
    for (double value : chain) {
        const double delta = value - chain_mean;
        var += delta * delta;
    }
    var /= static_cast<double>(n);
    if (var <= 0.0 || !std::isfinite(var)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double tau = 1.0;
    const std::size_t max_lag = n / 2;
    for (std::size_t lag = 1; lag <= max_lag; ++lag) {
        double cov = 0.0;
        for (std::size_t i = 0; i + lag < n; ++i) {
            cov += (chain[i] - chain_mean) * (chain[i + lag] - chain_mean);
        }
        cov /= static_cast<double>(n - lag);
        const double rho = cov / var;
        if (!std::isfinite(rho) || rho <= 0.0) {
            break;
        }
        tau += 2.0 * rho;
    }
    return tau;
}

double autocorr_time_for_parameter(const std::vector<std::vector<double>>& chains) {
    std::vector<double> taus;
    for (const auto& chain : chains) {
        const double tau = autocorr_time_for_chain(chain);
        if (std::isfinite(tau)) {
            taus.push_back(tau);
        }
    }
    return mean(taus);
}

void write_mcmc_diagnostics(
    const std::vector<MCMCSample>& samples,
    int n_steps,
    int n_walkers,
    int n_ensembles,
    int n_threads,
    double runtime_seconds,
    const std::string& path
) {
    constexpr double burnin_fraction = 0.25;
    const int burnin_step = static_cast<int>(burnin_fraction * static_cast<double>(n_steps));

    std::vector<int> accepted_by_ensemble(static_cast<std::size_t>(n_ensembles), 0);
    std::vector<int> count_by_ensemble(static_cast<std::size_t>(n_ensembles), 0);
    int accepted_total = 0;

    for (const auto& sample : samples) {
        if (sample.ensemble < 0 || sample.ensemble >= n_ensembles) {
            continue;
        }
        count_by_ensemble[static_cast<std::size_t>(sample.ensemble)] += 1;
        if (sample.accepted) {
            accepted_by_ensemble[static_cast<std::size_t>(sample.ensemble)] += 1;
            accepted_total += 1;
        }
    }

    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Nao foi possivel escrever: " + path);
    }

    out << std::setprecision(16);
    out << "# MCMC diagnostics\n";
    out << "n_steps " << n_steps << '\n';
    out << "n_walkers " << n_walkers << '\n';
    out << "n_ensembles " << n_ensembles << '\n';
    out << "n_threads " << n_threads << '\n';
    out << "runtime_seconds " << runtime_seconds << '\n';
    out << "runtime_minutes " << runtime_seconds / 60.0 << '\n';
    out << "burnin_fraction " << burnin_fraction << '\n';
    out << "burnin_steps " << burnin_step << '\n';
    out << "samples_total " << samples.size() << '\n';
    out << "samples_used_for_rhat_autocorr " << std::max(0, n_steps - burnin_step) * n_walkers * n_ensembles << '\n';

    std::vector<double> ensemble_acceptance_rates;
    ensemble_acceptance_rates.reserve(static_cast<std::size_t>(n_ensembles));
    for (int ensemble = 0; ensemble < n_ensembles; ++ensemble) {
        const int total = count_by_ensemble[static_cast<std::size_t>(ensemble)];
        const int accepted = accepted_by_ensemble[static_cast<std::size_t>(ensemble)];
        if (total > 0) {
            ensemble_acceptance_rates.push_back(static_cast<double>(accepted) / static_cast<double>(total));
        }
    }

    const auto minmax_acceptance = std::minmax_element(ensemble_acceptance_rates.begin(), ensemble_acceptance_rates.end());
    out << "acceptance_rate_global "
        << (samples.empty() ? std::numeric_limits<double>::quiet_NaN()
                            : static_cast<double>(accepted_total) / static_cast<double>(samples.size()))
        << '\n';
    out << "acceptance_rate_mean_ensembles " << mean(ensemble_acceptance_rates) << '\n';
    out << "acceptance_rate_min_ensemble "
        << (ensemble_acceptance_rates.empty() ? std::numeric_limits<double>::quiet_NaN() : *minmax_acceptance.first)
        << '\n';
    out << "acceptance_rate_max_ensemble "
        << (ensemble_acceptance_rates.empty() ? std::numeric_limits<double>::quiet_NaN() : *minmax_acceptance.second)
        << "\n\n";

    const std::array<const char*, 3> names = {"log10_lambda", "alpha_U", "w_U"};
    out << "# Gelman-Rubin R_hat/PSRF computed across independent ensembles after burn-in.\n";
    out << "# parameter gelman_rubin_R_hat autocorr_time_steps\n";
    for (int p = 0; p < 3; ++p) {
        std::vector<std::vector<double>> chains(static_cast<std::size_t>(n_ensembles));
        for (const auto& sample : samples) {
            if (sample.step < burnin_step || sample.ensemble < 0 || sample.ensemble >= n_ensembles) {
                continue;
            }
            const double value = sample_value(sample, p);
            if (std::isfinite(value)) {
                chains[static_cast<std::size_t>(sample.ensemble)].push_back(value);
            }
        }

        const std::size_t min_size = std::min_element(
            chains.begin(),
            chains.end(),
            [](const auto& a, const auto& b) { return a.size() < b.size(); }
        )->size();

        for (auto& chain : chains) {
            if (chain.size() > min_size) {
                chain.resize(min_size);
            }
        }

        out << names[static_cast<std::size_t>(p)] << ' '
            << rhat_for_parameter(chains) << ' '
            << autocorr_time_for_parameter(chains) << '\n';
    }
}

void write_likelihood_breakdown(
    const MCMCSample& best,
    const EOS& eos,
    const std::vector<double>& rho_c,
    const ObservationConfig& obs,
    const std::string& profile,
    const std::string& path
) {
    const std::array<double, 3> theta{best.log10_lambda, best.alpha_U, best.w_U};
    const LikelihoodBreakdown breakdown = log_likelihood_breakdown(theta, eos, rho_c, obs);

    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Nao foi possivel escrever: " + path);
    }

    out << std::setprecision(16);
    out << "# Likelihood breakdown for best post-burn-in posterior sample\n";
    out << "profile " << profile << '\n';
    out << "log10_lambda " << best.log10_lambda << '\n';
    out << "lambda " << best.lambda << '\n';
    out << "alpha_U " << best.alpha_U << '\n';
    out << "w_U " << best.w_U << '\n';
    out << "prior " << breakdown.prior << '\n';
    out << "logL_ligo_virgo_raw " << breakdown.ligo_virgo << '\n';
    out << "logL_nicer_j0740_raw " << breakdown.nicer_j0740 << '\n';
    out << "logL_nicer_j1231_raw " << breakdown.nicer_j1231 << '\n';
    out << "constraint_penalty " << breakdown.constraint_penalty << '\n';
    out << "weight_ligo_virgo " << obs.weight_ligo_virgo << '\n';
    out << "weight_nicer_j0740 " << obs.weight_nicer_j0740 << '\n';
    out << "weight_nicer_j1231 " << obs.weight_nicer_j1231 << '\n';
    out << "support_coverage_weight " << obs.support_coverage_weight << '\n';
    out << "logposterior_total " << breakdown.total << '\n';
    out << "Mmax " << breakdown.diag.Mmax << '\n';
    out << "Rmax " << breakdown.diag.Rmax << '\n';
    out << "R14 " << breakdown.diag.R14 << '\n';
    out << "Cmax " << breakdown.diag.Cmax << '\n';
    out << "Mmax_min " << obs.Mmax_min << '\n';
    out << "Mmax_sigma_below " << obs.Mmax_sigma_below << '\n';
    out << "Cmax_max " << obs.Cmax_max << '\n';
    out << "Cmax_sigma_above " << obs.Cmax_sigma_above << '\n';
    out << "radius_bandwidth_km " << obs.radius_bandwidth_km << '\n';
    out << "nicer_radius_bandwidth_km " << obs.nicer_radius_bandwidth_km << '\n';
    out << "prior_log10_lambda " << obs.log10_lambda_min << ' ' << obs.log10_lambda_max << '\n';
    out << "prior_log10_lambda_gaussian "
        << (obs.use_gaussian_log10_lambda_prior ? 1 : 0) << ' '
        << obs.log10_lambda_mu << ' '
        << obs.log10_lambda_sigma << '\n';
    out << "prior_alpha_U " << obs.alpha_U_min << ' ' << obs.alpha_U_max << '\n';
    out << "prior_alpha_U_gaussian "
        << (obs.use_gaussian_alpha_U_prior ? 1 : 0) << ' '
        << obs.alpha_U_mu << ' '
        << obs.alpha_U_sigma << '\n';
    out << "fixed_alpha_U "
        << (obs.fix_alpha_U ? 1 : 0) << ' '
        << obs.fixed_alpha_U << '\n';
    out << "prior_w_U " << obs.w_U_min << ' ' << obs.w_U_max << '\n';
    out << "prior_w_U_gaussian "
        << (obs.use_gaussian_w_U_prior ? 1 : 0) << ' '
        << obs.w_U_mu << ' '
        << obs.w_U_sigma << '\n';
    out << "fixed_w_U "
        << (obs.fix_w_U ? 1 : 0) << ' '
        << obs.fixed_w_U << '\n';
}

void write_mass_radius_posterior(
    const std::vector<MCMCSample>& samples,
    const EOS& eos,
    const std::vector<double>& rho_c,
    const std::string& path,
    int max_seq = 1000
) {
    constexpr double burnin_fraction = 0.25;
    const int max_sequences = max_seq;
    constexpr int n_mass = 80;

    const int max_step = samples.empty() ? 0 : samples.back().step;
    const int burnin_step = static_cast<int>(burnin_fraction * static_cast<double>(max_step + 1));

    std::vector<const MCMCSample*> posterior;
    for (const auto& sample : samples) {
        if (sample.step >= burnin_step && std::isfinite(sample.logposterior)) {
            posterior.push_back(&sample);
        }
    }

    std::vector<std::size_t> indices(posterior.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::mt19937 rng(42);
    std::shuffle(indices.begin(), indices.end(), rng);
    const std::size_t n_draw = std::min<std::size_t>(static_cast<std::size_t>(max_sequences), indices.size());
    indices.resize(n_draw);

    std::vector<std::array<double, n_mass>> radii_by_sequence;
    radii_by_sequence.reserve(n_draw);

    const double m_min = 0.8;
    const double m_max = 2.5;
    std::array<double, n_mass> mass_grid{};
    for (int i = 0; i < n_mass; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(n_mass - 1);
        mass_grid[static_cast<std::size_t>(i)] = m_min + t * (m_max - m_min);
    }

    for (const std::size_t idx : indices) {
        const auto& sample = *posterior[idx];
        const BraneParams params{sample.lambda, sample.alpha_U, sample.w_U};
        const Sequence seq = mass_radius_sequence(rho_c, eos, params, "posterior", 80.0, 0.08);

        std::array<double, n_mass> radii{};
        for (int j = 0; j < n_mass; ++j) {
            radii[static_cast<std::size_t>(j)] = radius_at_mass(seq, mass_grid[static_cast<std::size_t>(j)]);
        }
        radii_by_sequence.push_back(radii);
    }

    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Nao foi possivel escrever: " + path);
    }

    out << "# M_sun R_q025 R_q16 R_median R_mean R_q84 R_q975 n\n";
    out << std::setprecision(16);
    for (int j = 0; j < n_mass; ++j) {
        std::vector<double> radii;
        for (const auto& row : radii_by_sequence) {
            if (std::isfinite(row[static_cast<std::size_t>(j)])) {
                radii.push_back(row[static_cast<std::size_t>(j)]);
            }
        }
        const double coverage = radii_by_sequence.empty()
            ? 0.0
            : static_cast<double>(radii.size()) / static_cast<double>(radii_by_sequence.size());
        if (radii.size() < 10 || coverage < 0.5) {
            continue;
        }
        out << mass_grid[static_cast<std::size_t>(j)] << ' '
            << quantile(radii, 0.025) << ' '
            << quantile(radii, 0.16) << ' '
            << quantile(radii, 0.50) << ' '
            << mean(radii) << ' '
            << quantile(radii, 0.84) << ' '
            << quantile(radii, 0.975) << ' '
            << radii.size() << '\n';
    }
}

bool sufficiently_different(const MCMCSample& a, const MCMCSample& b) {
    return std::abs(a.log10_lambda - b.log10_lambda) >= 0.05
        || std::abs(a.alpha_U - b.alpha_U) >= 0.02
        || std::abs(a.w_U - b.w_U) >= 0.03;
}

std::string best_curve_label(const MCMCSample& sample, int rank) {
    std::ostringstream label;
    label << "best " << rank
          << ": lambda=" << std::scientific << std::setprecision(2) << sample.lambda
          << ", alpha_U=" << std::fixed << std::setprecision(3) << sample.alpha_U
          << ", w_U=" << std::fixed << std::setprecision(3) << sample.w_U;
    return label.str();
}

void write_best_stable_curves(
    const std::vector<MCMCSample>& samples,
    const EOS& eos,
    const std::vector<double>& rho_c,
    const std::string& output_dir
) {
    constexpr double burnin_fraction = 0.25;
    const int max_step = samples.empty() ? 0 : samples.back().step;
    const int burnin_step = static_cast<int>(burnin_fraction * static_cast<double>(max_step + 1));

    std::vector<const MCMCSample*> ranked;
    for (const auto& sample : samples) {
        if (sample.step >= burnin_step && std::isfinite(sample.logposterior)) {
            ranked.push_back(&sample);
        }
    }

    std::sort(ranked.begin(), ranked.end(), [](const auto* a, const auto* b) {
        return a->logposterior > b->logposterior;
    });

    std::vector<const MCMCSample*> selected;
    for (const auto* sample : ranked) {
        bool different = true;
        for (const auto* previous : selected) {
            if (!sufficiently_different(*sample, *previous)) {
                different = false;
                break;
            }
        }
        if (different) {
            selected.push_back(sample);
        }
        if (selected.size() == 2) {
            break;
        }
    }

    for (int rank = 1; rank <= 2; ++rank) {
        const std::string path = (std::filesystem::path(output_dir) / ("best_curve_" + std::to_string(rank) + ".dat")).string();
        if (static_cast<std::size_t>(rank) > selected.size()) {
            std::ofstream out(path);
            out << "# label best " << rank << " unavailable\n";
            out << "# rho_c_cgs rho_c_geom p_c_geom R_km M_km M_sun compactness lambda alpha_U w_U\n";
            continue;
        }

        const auto& sample = *selected[static_cast<std::size_t>(rank - 1)];
        const BraneParams params{sample.lambda, sample.alpha_U, sample.w_U};
        const Sequence seq = mass_radius_sequence(rho_c, eos, params, best_curve_label(sample, rank), 80.0, 0.05);
        write_stable_sequence_dat(seq, path);
    }
}

std::vector<MCMCSample> read_mcmc_dat(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Nao foi possivel ler: " + path);
    }

    std::vector<MCMCSample> samples;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }

        MCMCSample sample;
        int accepted = 0;
        std::istringstream iss(line);
        iss >> sample.ensemble
            >> sample.step
            >> sample.walker
            >> sample.log10_lambda
            >> sample.lambda
            >> sample.alpha_U
            >> sample.w_U
            >> sample.logposterior
            >> sample.diag.Mmax
            >> sample.diag.Rmax
            >> sample.diag.R14
            >> sample.diag.Cmax
            >> sample.diag.rhoc_Mmax
            >> accepted;
        if (!iss) {
            continue;
        }
        sample.diag.valid = std::isfinite(sample.diag.Mmax);
        sample.accepted = accepted != 0;
        samples.push_back(sample);
    }

    return samples;
}

const MCMCSample* best_post_burnin_sample(const std::vector<MCMCSample>& samples, int n_steps) {
    constexpr double burnin_fraction = 0.25;
    const int burnin_step = static_cast<int>(burnin_fraction * static_cast<double>(n_steps));

    const MCMCSample* best = nullptr;
    for (const auto& sample : samples) {
        if (sample.step < burnin_step || !std::isfinite(sample.logposterior)) {
            continue;
        }
        if (best == nullptr || sample.logposterior > best->logposterior) {
            best = &sample;
        }
    }
    return best;
}

void run_breakdown_from_chain(
    const EOS& eos,
    const std::vector<double>& rho_c,
    const std::string& profile,
    int n_steps,
    const std::string& output_dir
) {
    const ObservationConfig obs = observation_profile(profile);
    const std::string chain_path = (std::filesystem::path(output_dir) / "mcmc_chain.dat").string();
    const auto samples = read_mcmc_dat(chain_path);
    const MCMCSample* best = best_post_burnin_sample(samples, n_steps);
    if (best == nullptr) {
        throw std::runtime_error("Nenhuma amostra pos-burn-in valida encontrada em " + chain_path);
    }

    write_likelihood_breakdown(*best, eos, rho_c, obs, profile, (std::filesystem::path(output_dir) / "likelihood_breakdown.dat").string());
    std::cout << "Likelihood breakdown atualizado para o melhor ponto pos-burn-in:\n"
              << "  log10(lambda) = " << best->log10_lambda << '\n'
              << "  lambda        = " << best->lambda << '\n'
              << "  alpha_U       = " << best->alpha_U << '\n'
              << "  w_U           = " << best->w_U << '\n'
              << "  logpost       = " << best->logposterior << '\n'
              << "  Mmax          = " << best->diag.Mmax << " M_sun\n"
              << "  R1.4          = " << best->diag.R14 << " km\n";
}

void run_chain(
    const EOS& eos,
    const std::vector<double>& rho_c,
    int n_steps,
    int n_walkers,
    int n_ensembles,
    int n_threads,
    const std::string& profile,
    const std::string& output_dir
) {
    std::filesystem::create_directories(output_dir);

    MCMCConfig config;
    config.n_steps = n_steps;
    config.n_walkers = n_walkers;
    config.n_ensembles = n_ensembles;
    config.n_threads = n_threads;
    const ObservationConfig obs = observation_profile(profile);
    // Initialization now comes from the profile (obs.init_*), never from a
    // constant hard-coded here: a run must be fully determined by
    // (profile, seed, n_steps, n_walkers, n_ensembles). See docs/PHASE0.md.
    if (obs.override_init) {
        config.theta0 = obs.init_center;
        config.initial_width = obs.init_width;
    }

    write_run_provenance(profile, obs, config, eos.name(), rho_c, output_dir);

    std::cout << "Rodando " << n_ensembles << " ensembles independentes tipo emcee com "
              << n_walkers << " walkers, " << n_steps << " passos e "
              << n_threads << " threads, profile=" << profile << "...\n";
    const auto start = std::chrono::steady_clock::now();
    const auto samples = run_mcmc(eos, rho_c, obs, config);
    const auto end = std::chrono::steady_clock::now();
    const double runtime_seconds = std::chrono::duration<double>(end - start).count();
    std::cout << "Escrevendo cadeia e diagnosticos...\n";
    write_mcmc_dat(samples, (std::filesystem::path(output_dir) / "mcmc_chain.dat").string());
    write_mcmc_diagnostics(samples, n_steps, n_walkers, n_ensembles, n_threads, runtime_seconds, (std::filesystem::path(output_dir) / "mcmc_diagnostics.dat").string());
    std::cout << "Gerando bandas massa-raio posteriores...\n";
    write_mass_radius_posterior(samples, eos, rho_c, (std::filesystem::path(output_dir) / "mass_radius_posterior.dat").string());
    std::cout << "Gerando curvas auxiliares...\n";
    write_best_stable_curves(samples, eos, rho_c, output_dir);
    const Sequence gr_seq = mass_radius_sequence(
        rho_c,
        eos,
        {std::numeric_limits<double>::infinity(), 0.0, 0.0},
        "GR",
        80.0,
        0.05
    );
    write_stable_sequence_dat(gr_seq, (std::filesystem::path(output_dir) / "sequence_gr.dat").string());

    const int accepted = static_cast<int>(std::count_if(samples.begin(), samples.end(), [](const MCMCSample& s) {
        return s.accepted;
    }));
    const MCMCSample* best = best_post_burnin_sample(samples, n_steps);

    std::cout << "Taxa de aceitacao = " << static_cast<double>(accepted) / samples.size() << '\n';
    std::cout << "Tempo total do MCMC = " << runtime_seconds << " s (" << runtime_seconds / 60.0 << " min)\n";
    std::cout << "Resultados em " << output_dir << '\n';
    if (best != nullptr) {
        write_likelihood_breakdown(*best, eos, rho_c, obs, profile, (std::filesystem::path(output_dir) / "likelihood_breakdown.dat").string());
        std::cout << "Melhor ponto pos-burn-in:\n"
                  << "  log10(lambda) = " << best->log10_lambda << '\n'
                  << "  lambda        = " << best->lambda << '\n'
                  << "  alpha_U       = " << best->alpha_U << '\n'
                  << "  w_U           = " << best->w_U << '\n'
                  << "  logpost       = " << best->logposterior << '\n'
                  << "  Mmax          = " << best->diag.Mmax << " M_sun\n"
                  << "  R1.4          = " << best->diag.R14 << " km\n";
    }
    // Last, so the figure covers the posterior bands and auxiliary curves too.
    append_peak_rss(output_dir, runtime_seconds, samples.size(), 1000);
}

} // namespace

int main(int argc, char** argv) {
    const std::string mode = argc > 1 ? argv[1] : "all";

    if (mode == "help" || mode == "--help" || mode == "-h") {
        print_usage(argv[0]);
        return 0;
    }

    try {
        std::string eos_path = "data/SLy.txt";
        std::string output_dir = "output/data";
        if (mode == "sequence") {
            eos_path = argc > 2 ? argv[2] : eos_path;
            output_dir = argc > 3 ? argv[3] : output_dir;
        } else if (mode == "mcmc" || mode == "all") {
            eos_path = argc > 7 ? argv[7] : eos_path;
            output_dir = argc > 8 ? argv[8] : output_dir;
        } else if (mode == "breakdown") {
            eos_path = argc > 4 ? argv[4] : eos_path;
            output_dir = argc > 5 ? argv[5] : output_dir;
        } else if (mode == "reband") {
            eos_path = argc > 4 ? argv[4] : eos_path;
            output_dir = argc > 5 ? argv[5] : output_dir;
        } else if (mode == "sequence-custom") {
            eos_path = argc > 6 ? argv[6] : eos_path;
        } else if (mode == "profile") {
            eos_path = argc > 7 ? argv[7] : eos_path;
        } else if (mode == "serve") {
            eos_path = argc > 3 ? argv[3] : eos_path;
        } else if (mode == "eosinfo") {
            eos_path = argc > 2 ? argv[2] : eos_path;
        }

        // Integration configuration. Defaults are the values validated in
        // docs/PHASE_C1.md; the env overrides exist for sensitivity sweeps, and the
        // RESOLVED values are written to effective_config.dat whatever set them, so
        // provenance does not depend on how the job was launched.
        IntegratorConfig icfg;
        if (const char* v = std::getenv("BRANA_RTOL")) { icfg.rtol = std::atof(v); }
        if (const char* v = std::getenv("BRANA_ATOL")) { icfg.atol = std::atof(v); }
        if (const char* v = std::getenv("BRANA_SIGMA_SURFACE")) { icfg.sigma_surface = std::atof(v); }
        if (const char* v = std::getenv("BRANA_P_STOP_FRAC")) { icfg.p_stop_frac = std::atof(v); }
        if (const char* v = std::getenv("BRANA_RHO_JOIN_CGS")) { icfg.rho_join_cgs = std::atof(v); }
        if (std::getenv("BRANA_NO_LIMITER") != nullptr) { icfg.use_surface_limiter = false; }
        if (std::getenv("BRANA_TIDAL_IN_ERRCTL") != nullptr) { icfg.tidal_in_error_control = true; }
        set_integrator_config(icfg);
        const auto eos = TabulatedEOS::from_file_cgs(eos_path, icfg.rho_join_cgs);
        const int n_rho = mode == "sequence-custom" && argc > 7 ? std::stoi(argv[7]) : 120;
        const double custom_dr = mode == "sequence-custom" && argc > 8 ? std::stod(argv[8]) : 0.05;
        if (n_rho < 20) {
            throw std::runtime_error("n_rho deve ser pelo menos 20.");
        }
        if (!(custom_dr > 0.0)) {
            throw std::runtime_error("dr deve ser positivo.");
        }
        const auto rho_c = central_densities_tabulated(eos, n_rho);

        if (mode != "serve" && mode != "eosinfo") {   // stdout belongs to the protocol
            std::cout << eos.name() << '\n';
            std::cout << "output_dir=" << output_dir << '\n';
            std::cout << "rho_c usado de " << rho_geom_km2_to_cgs(rho_c.front()) << " ate "
                      << rho_geom_km2_to_cgs(rho_c.back()) << " g/cm^3\n";
        }

        if (mode == "sequence" || mode == "all") {
            run_sequences(eos, rho_c, output_dir);
        }
        if (mode == "sequence-custom") {
            if (argc < 5) {
                print_usage(argv[0]);
                return 1;
            }
            const double lambda_brane = std::stod(argv[2]);
            const double alpha_U = std::stod(argv[3]);
            const double w_U = std::stod(argv[4]);
            const std::string path = argc > 5 ? argv[5] : "output/data/sequence_modified_custom.dat";
            const auto closure = argc > 9
                ? static_cast<WeylClosure>(std::stoi(argv[9])) : WeylClosure::Original;
            const double rho0_cgs = argc > 10 ? std::stod(argv[10]) : 0.0;
            const double n_sup = argc > 11 ? std::stod(argv[11]) : 3.0;
            run_custom_sequence(eos, rho_c, lambda_brane, alpha_U, w_U, path, custom_dr,
                                closure, rho0_cgs, n_sup);
        }
        if (mode == "profile") {
            if (argc < 7) {
                print_usage(argv[0]);
                return 1;
            }
            const double lambda_brane = std::stod(argv[2]);
            const double alpha_U = std::stod(argv[3]);
            const double w_U = std::stod(argv[4]);
            const double rho_c_cgs = std::stod(argv[5]);
            const std::string path = argv[6];
            const double prof_dr = argc > 8 ? std::stod(argv[8]) : 0.01;
            const auto closure = argc > 9
                ? static_cast<WeylClosure>(std::stoi(argv[9])) : WeylClosure::Original;
            const double rho0_cgs = argc > 10 ? std::stod(argv[10]) : 0.0;
            const double n_sup = argc > 11 ? std::stod(argv[11]) : 3.0;
            run_profile(eos, lambda_brane, alpha_U, w_U, rho_c_cgs, path, prof_dr,
                        closure, rho0_cgs, n_sup);
        }
        if (mode == "mcmc" || mode == "all") {
            const int n_steps = argc > 2 ? std::stoi(argv[2]) : 800;
            const int n_walkers = argc > 3 ? std::stoi(argv[3]) : 32;
            const int n_ensembles = argc > 4 ? std::stoi(argv[4]) : 8;
            const int n_threads = argc > 5 ? std::stoi(argv[5]) : 8;
            const std::string profile = argc > 6 ? argv[6] : "default";
            run_chain(eos, rho_c, n_steps, n_walkers, n_ensembles, n_threads, profile, output_dir);
        }
        if (mode == "breakdown") {
            const std::string profile = argc > 2 ? argv[2] : "perturbative-gr";
            const int n_steps = argc > 3 ? std::stoi(argv[3]) : 3000;
            run_breakdown_from_chain(eos, rho_c, profile, n_steps, output_dir);
        }
        if (mode == "reband") {
            const std::string chain_path = argc > 2 ? argv[2] : "output/data/mcmc_chain.dat";
            const int max_seq = argc > 3 ? std::stoi(argv[3]) : 1000;
            std::cout << "Lendo cadeia MCMC de " << chain_path << "...\n";
            const auto samples = read_mcmc_dat(chain_path);
            std::cout << "Cadeia carregada: " << samples.size() << " amostras.\n";
            std::cout << "Gerando banda M-R com " << max_seq << " sequencias...\n";
            std::filesystem::create_directories(output_dir);
            const std::string band_path = (std::filesystem::path(output_dir) / "mass_radius_posterior.dat").string();
            write_mass_radius_posterior(samples, eos, rho_c, band_path, max_seq);
            std::cout << "Banda M-R em " << band_path << '\n';
        }
        if (mode == "eosinfo") {
            // Atmosphere diagnostic: parameters and continuity across the join.
            const double rj = eos.atmosphere_rho_join();
            std::cout << std::setprecision(10);
            std::cout << "atmosphere " << (eos.has_atmosphere() ? 1 : 0) << '\n';
            if (eos.has_atmosphere()) {
                std::cout << "rho_join_cgs " << rho_geom_km2_to_cgs(rj) << '\n';
                std::cout << "gamma_atm " << eos.atmosphere_gamma() << '\n';
                std::cout << "table_rho_min_cgs " << rho_geom_km2_to_cgs(eos.rho_max()) << ' '
                          << "(topo)\n";
                for (double f : {0.99, 0.999, 1.0, 1.001, 1.01}) {
                    const double r = rj * f;
                    std::cout << "probe " << f << ' ' << eos.p_of_rho(r) << ' '
                              << eos.cs2_of_rho(r) << '\n';
                }
            }
            return 0;
        }
        if (mode == "serve") {
            // args: serve <profile> <eos.txt>. Quiet on stdout: the protocol owns it.
            const std::string profile = argc > 2 ? argv[2] : "perturbative-gr";
            run_serve(eos, rho_c, observation_profile(profile));
            return 0;
        }
        if (mode != "sequence" && mode != "sequence-custom" && mode != "mcmc" && mode != "breakdown" && mode != "reband" && mode != "profile" && mode != "serve" && mode != "eosinfo" && mode != "all") {
            print_usage(argv[0]);
            return 1;
        }
    } catch (const std::exception& exc) {
        std::cerr << "Erro: " << exc.what() << '\n';
        return 1;
    }

    return 0;
}
