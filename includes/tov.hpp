#pragma once

#include "eos.hpp"

#include <limits>
#include <string>
#include <vector>

namespace brana {

// Which closure is used for the Weyl sector. The selection criterion is a stated
// principle, not the score on a table. See docs/PHASE_C1.md.
//
// SELF-SIMILARITY. The projected Weyl tensor decomposes into a nonlocal energy
// density U and a nonlocal anisotropic stress P (Germani & Maartens 2001), and
// the quantities this code integrates are
//     rho_U = 6 U / (kappa^4 lambda),   p_U = (2 U + 4 P) / (kappa^4 lambda),
// so that, exactly,
//     w_U == p_U / rho_U == 1/3 + (2/3) (P / U).
// Therefore w_U is not a free phenomenological knob: it MEASURES the dark
// anisotropic stress, and w_U = 1/3 is the isotropic (traceless, "dark
// radiation") sector. For w_U to have one geometric meaning throughout the star,
// P/U must not depend on radius. We require that.
//
// Consequence, and the reason V1 is the production closure:
//   Original : rho_U = alpha_U rho,        p_U = w_U rho_U
//              P/U constant -> self-similar. But rho_U(R) = alpha_U rho_surf != 0
//              on a tabulated crust, so p_r,tot(R) != 0 and the junction residual
//              DIVERGES under refinement in dr (see docs/PHASE_C1.md).
//   V1       : rho_U = alpha_U rho f(rho), p_U = w_U rho_U,  f = 1/(1+(rho0/rho)^n)
//              Any suppression of this form is automatically self-similar: f
//              multiplies U and P alike, so P/U is untouched and w_U keeps its
//              meaning in the crust as well as the core. It also drives rho_U(R)
//              to ~1e-23 of its unsuppressed value, which closes the junction.
//   V2       : rho_U = alpha_U rho,        p_U = w_U p
//              FAILS self-similarity: p_U/rho_U = w_U p/(alpha_U rho) varies with
//              density, so P/U becomes radius-dependent and w_U loses its
//              geometric referent. Kept runnable for the record, not for
//              production.
// The V1 suppression parameters are FIXED (rho0 = 2 rho_sat, n = 3), never sampled.
enum class WeylClosure { Original = 0, V1Suppressed = 1, V2PressureTracking = 2 };

struct BraneParams {
    double lambda_brane = std::numeric_limits<double>::infinity();
    double alpha_U = 0.0;
    double w_U = 0.0;
    WeylClosure closure = WeylClosure::Original;
    double rho0_geom = 0.0;   // V1 suppression scale (geometric units); 0 disables
    double n_suppress = 3.0;  // V1 suppression exponent
    // Phase 2 bracket, parametrised: delta U / U = beta * delta rho / rho.
    //   beta = 0 -> rigid Weyl   (dp/deps = cs2, matter alone responds)
    //   beta = 1 -> comoving     (dp/deps = cs2_eff, the full effective EoS)
    // In between, dp/deps = (cs2 + beta w alpha g)/(1 + beta alpha g), which is
    // exactly cs2_eff evaluated with alpha -> beta*alpha. The 1/lambda terms are not
    // scaled by beta; at the large tensions sampled here they are negligible.
    double tidal_beta = 0.0;
};

// Weyl energy density for the active closure.
double rho_weyl(double rho, const BraneParams& params);

struct Star {
    double y_R = std::numeric_limits<double>::quiet_NaN();   // r H'/H at the surface
    double R = std::numeric_limits<double>::quiet_NaN();
    double M_km = std::numeric_limits<double>::quiet_NaN();
    double M_sun = std::numeric_limits<double>::quiet_NaN();
    double compactness = std::numeric_limits<double>::quiet_NaN();
    double rho_c = std::numeric_limits<double>::quiet_NaN();
    double rho_c_cgs = std::numeric_limits<double>::quiet_NaN();
    double p_c = std::numeric_limits<double>::quiet_NaN();
    BraneParams params;
};

// One accepted RK step of the interior solution.
struct ProfilePoint {
    double r = 0.0;
    double m = 0.0;         // enclosed mass in km (uses rho_tot, as the TOV does)
    double p = 0.0;         // matter pressure, geometric
    double rho = 0.0;       // matter energy density, geometric
    double rho_tot = 0.0;   // effective density entering dm/dr
    double p_r_tot = 0.0;   // effective radial pressure entering dp/dr
    double cs2 = 0.0;       // matter dp/drho
    double cs2_eff = 0.0;   // d p_r_tot / d rho_tot along the sequence
};

struct Sequence {
    std::string label;
    std::vector<Star> stars;
};

// Numerical configuration of the stellar integration. These were constexpr and
// environment variables during Phase C7 development; a result that depends on how
// the job was launched is exactly the Phase 0 problem, so the resolved values are
// recorded in effective_config.dat whatever set them.
//
// Defaults are the values validated in docs/PHASE_C1.md: R_1.4 independent of dr to
// 3e-8 km across dr = 0.08 -> 0.005, against 0.068 km before.
struct IntegratorConfig {
    double rtol = 1.0e-7;
    double atol = 1.0e-20;        // was 1e-12, a floor larger than surface pressures
    double sigma_surface = 0.2;   // a step may consume at most this fraction of p
    double p_stop_frac = 1.0e-16; // surface floor relative to p_c, needed because the
                                  // limiter makes p decay geometrically toward zero
    double rho_join_cgs = 0.0;    // 0 = automatic atmosphere rule; <0 disables it
    bool use_surface_limiter = true;
    // Including y in the step control changes the accepted steps and therefore every
    // archived chain, so it is off by default and exists to measure whether the
    // structure steps are fine enough for the perturbation.
    bool tidal_in_error_control = false;
};

const IntegratorConfig& integrator_config();
void set_integrator_config(const IntegratorConfig& cfg);

struct Diagnostics {
    bool valid = false;
    double Mmax = std::numeric_limits<double>::quiet_NaN();
    double Rmax = std::numeric_limits<double>::quiet_NaN();
    double rhoc_Mmax = std::numeric_limits<double>::quiet_NaN();
    double Cmax = std::numeric_limits<double>::quiet_NaN();
    double R14 = std::numeric_limits<double>::quiet_NaN();
};

double rho_eff_weyl(double rho, const BraneParams& params);
double p_eff_weyl(double p, double rho, const BraneParams& params);

Star integrate_star(
    double rho_c,
    const EOS& eos,
    const BraneParams& params,
    double r_max = 80.0,
    double dr = 0.05,
    double r0 = 1.0e-5,
    std::vector<ProfilePoint>* profile = nullptr
);

// d p_r,tot / d rho_tot for the effective source, at fixed (lambda, alpha, w).
// In the large-tension limit this reduces to (cs2 + w*alpha) / (1 + alpha).
double cs2_eff_weyl(double rho, double p, double cs2, const BraneParams& params);

Sequence mass_radius_sequence(
    const std::vector<double>& rho_c_values,
    const EOS& eos,
    const BraneParams& params,
    const std::string& label,
    double r_max = 80.0,
    double dr = 0.05
);

Diagnostics sequence_diagnostics(const Sequence& seq);
std::vector<double> logspace(double log10_min, double log10_max, int n);

void write_sequence_dat(const Sequence& seq, const std::string& path);
void write_stable_sequence_dat(const Sequence& seq, const std::string& path);

} // namespace brana
