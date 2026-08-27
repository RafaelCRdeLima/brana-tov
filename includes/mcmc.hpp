#pragma once

#include "eos.hpp"
#include "tov.hpp"

#include <array>
#include <random>
#include <string>
#include <vector>

namespace brana {

struct ObservationConfig {
    std::string ligo_virgo_samples_path = "data/ligo_virgo/Parametrized-EoS_maxmass_posterior_samples.dat";
    double radius_bandwidth_km = 0.70;
    std::string nicer_j0740_samples_path = "data/nicer/J0740/J0740_NICERXMM_full_mr.txt";
    std::string nicer_j1231_samples_path = "data/nicer/J1231/mr_samples_and_contours/PDTU_H_R1014/J1231_R1014_wmrsamples.txt";
    double nicer_radius_bandwidth_km = 0.70;
    int max_nicer_samples_per_source = 20000;
    bool use_ligo_virgo = true;
    bool use_nicer_j0740 = true;
    bool use_nicer_j1231 = true;
    double weight_ligo_virgo = 1.0;
    double weight_nicer_j0740 = 1.0;
    double weight_nicer_j1231 = 1.0;
    double support_coverage_weight = 1.0;
    // Floor under the coverage fraction before taking its log, an undocumented
    // constant chosen only to avoid log(0). Measured, it does NOT bind on the
    // path between the SLy modes: profiling log-posterior along that segment gives
    // the same barrier (16.41) at 1e-12 and at 1e-6, because coverage in the
    // valley falls only to ~25% of its modal value, nowhere near the floor. Kept
    // configurable so the knob is visible and its irrelevance is checkable rather
    // than assumed. See docs/PHASE_C8.md.
    double support_coverage_floor = 1.0e-12;
    double log10_lambda_min = 0.0;
    double log10_lambda_max = 6.0;
    bool use_gaussian_log10_lambda_prior = true;
    double log10_lambda_mu = 4.0;
    double log10_lambda_sigma = 1.5;
    double alpha_U_min = -2.0;
    double alpha_U_max = 0.5;
    bool use_gaussian_alpha_U_prior = false;
    double alpha_U_mu = 0.0;
    double alpha_U_sigma = 0.5;
    bool fix_alpha_U = false;
    double fixed_alpha_U = 0.0;
    double w_U_min = -2.0;
    double w_U_max = 1.0;
    bool use_gaussian_w_U_prior = false;
    double w_U_mu = 0.0;
    double w_U_sigma = 0.5;
    bool fix_w_U = false;
    double fixed_w_U = 0.0;
    double Mmax_min = 1.97;
    double Mmax_sigma_below = 0.05;
    double Cmax_max = 0.34;
    double Cmax_sigma_above = 0.03;

    // Walker initialization. Part of the profile so that a run is fully
    // determined by (profile, seed, n_steps, n_walkers, n_ensembles) and never
    // by which binary happened to be on disk. See docs/PHASE0.md.
    bool override_init = false;                       // if false, MCMCConfig defaults are kept
    std::array<double, 3> init_center = {3.0, -0.25, -0.1};
    std::array<double, 3> init_width = {0.35, 0.45, 0.45};

    // Weyl closure (Phase C1). Part of the profile for the same reason the
    // initialization is: the closure changes the TOV solution and therefore the
    // likelihood, so a V2 posterior is NOT the original posterior filtered.
    WeylClosure closure = WeylClosure::Original;
    double closure_rho0_cgs = 0.0;      // V1 suppression scale; 0 disables
    double closure_n_suppress = 3.0;    // V1 suppression exponent

    // Ablation switches (Phase 1). Defaults reproduce the production likelihood.
    bool use_support_cut = true;        // hard -inf if <10% of samples lie in the model mass support
    double support_cut_fraction = 0.1;
    bool use_constraint_penalties = true;  // Mmax_min / Cmax_max soft penalties
    // Hard rejection where the effective EoS is non-monotonic, d p_tot/d rho_tot < 0.
    // Hard and not a soft penalty on purpose: a non-monotonic p_tot(rho_tot) does not
    // define a constitutive relation at all, so there is no meaningful gradation of
    // "how little monotonic", and a soft version would reintroduce the ad hoc knob
    // this work criticises. The discarded prior fraction is a property of the EoS
    // microphysics -- it tracks min c_s^2 near rho_sat -- not of the closure.
    bool use_monotonicity_cut = false;
};

struct MCMCConfig {
    int n_steps = 800;
    int n_walkers = 32;
    int n_ensembles = 8;
    int n_threads = 8;
    double stretch_a = 2.0;
    std::array<double, 3> theta0 = {3.0, -0.25, -0.1};
    std::array<double, 3> initial_width = {0.35, 0.45, 0.45};
    unsigned int seed = 12345;
    bool show_progress = true;
};

struct MCMCSample {
    int ensemble = 0;
    int step = 0;
    int walker = 0;
    double log10_lambda = 0.0;
    double lambda = 0.0;
    double alpha_U = 0.0;
    double w_U = 0.0;
    double logposterior = 0.0;
    Diagnostics diag;
    bool accepted = false;
};

struct LikelihoodBreakdown {
    double prior = 0.0;
    double ligo_virgo = 0.0;
    double nicer_j0740 = 0.0;
    double nicer_j1231 = 0.0;
    double constraint_penalty = 0.0;
    double total = 0.0;
    Diagnostics diag;
};

double log_prior(const std::array<double, 3>& theta, const ObservationConfig& obs);
double log_likelihood(
    const std::array<double, 3>& theta,
    const EOS& eos,
    const std::vector<double>& rho_c_values,
    const ObservationConfig& obs,
    Diagnostics* diag_out = nullptr
);

LikelihoodBreakdown log_likelihood_breakdown(
    const std::array<double, 3>& theta,
    const EOS& eos,
    const std::vector<double>& rho_c_values,
    const ObservationConfig& obs
);

double log_posterior(
    const std::array<double, 3>& theta,
    const EOS& eos,
    const std::vector<double>& rho_c_values,
    const ObservationConfig& obs,
    Diagnostics* diag_out = nullptr
);

std::vector<MCMCSample> run_mcmc(
    const EOS& eos,
    const std::vector<double>& rho_c_values,
    const ObservationConfig& obs,
    const MCMCConfig& config
);

void write_mcmc_dat(const std::vector<MCMCSample>& samples, const std::string& path);

} // namespace brana
