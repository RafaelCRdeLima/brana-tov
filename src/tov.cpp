#include "tov.hpp"

#include "constants.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace brana {

namespace {
IntegratorConfig g_integrator_config{};
}

const IntegratorConfig& integrator_config() { return g_integrator_config; }
void set_integrator_config(const IntegratorConfig& cfg) { g_integrator_config = cfg; }

namespace {

struct State {
    double m = 0.0;
    double p = 0.0;
    // l = 2 metric perturbation, y = r H'/H, integrated alongside the structure so
    // that no interpolation of the background is needed (Hinderer 2008). It does not
    // feed back into (m, p) and, by default, does not enter the step-size control --
    // otherwise every archived chain would shift. See docs/PHASE_C7.md.
    double y = 2.0;
};

State operator+(const State& a, const State& b) {
    return {a.m + b.m, a.p + b.p, a.y + b.y};
}

State operator-(const State& a, const State& b) {
    return {a.m - b.m, a.p - b.p, a.y - b.y};
}

State operator*(double factor, const State& y) {
    return {factor * y.m, factor * y.p, factor * y.y};
}

State rhs(double r, const State& y, const EOS& eos, const BraneParams& params) {
    if (y.p <= 0.0) {
        return {};
    }

    const double rho = eos.rho_of_p(y.p);
    if (!std::isfinite(rho)) {
        return {};
    }

    const double rho_tot = rho_eff_weyl(rho, params);
    const double pr_tot = p_eff_weyl(y.p, rho, params);
    const double dm_dr = 4.0 * pi * r * r * rho_tot;
    const double denom = r * (r - 2.0 * y.m);
    double dp_dr = 0.0;

    if (denom > 0.0) {
        dp_dr = -(rho + y.p) * (y.m + 4.0 * pi * r * r * r * pr_tot) / denom;
    }

    // dy/dr from  r y' + y^2 + y F + r^2 Q = 0.  eps and p here are the TOTAL
    // effective source, because that is what the Einstein equations see; the closure
    // ambiguity enters only through dp/deps, selected by params.tidal_comoving.
    double dy_dr = 0.0;
    if (denom > 0.0 && r > 0.0) {
        const double cs2 = eos.cs2_of_rho(rho);
        double dpde = cs2;
        if (params.tidal_beta != 0.0) {
            BraneParams scaled = params;
            scaled.alpha_U *= params.tidal_beta;
            dpde = cs2_eff_weyl(rho, y.p, cs2, scaled);
        }
        const double one_minus = 1.0 - 2.0 * y.m / r;
        if (one_minus > 0.0 && std::isfinite(dpde) && dpde > 0.0) {
            const double elam = 1.0 / one_minus;
            const double F = elam * (1.0 - 4.0 * pi * r * r * (rho_tot - pr_tot));
            const double nup = 2.0 * (y.m + 4.0 * pi * r * r * r * pr_tot) / denom;
            const double r2Q =
                4.0 * pi * r * r * (5.0 * rho_tot + 9.0 * pr_tot + (rho_tot + pr_tot) / dpde) * elam
                - 6.0 * elam - (nup * r) * (nup * r);
            dy_dr = (-y.y * y.y - y.y * F - r2Q) / r;
        }
    }

    return {dm_dr, dp_dr, dy_dr};
}

struct RKStep {
    State y5;
    State error;
};

RKStep rk45_step(double r, const State& y, double h, const EOS& eos, const BraneParams& params) {
    const State k1 = rhs(r, y, eos, params);
    const State k2 = rhs(r + h * (1.0 / 5.0), y + h * ((1.0 / 5.0) * k1), eos, params);
    const State k3 = rhs(
        r + h * (3.0 / 10.0),
        y + h * ((3.0 / 40.0) * k1 + (9.0 / 40.0) * k2),
        eos,
        params
    );
    const State k4 = rhs(
        r + h * (4.0 / 5.0),
        y + h * ((44.0 / 45.0) * k1 + (-56.0 / 15.0) * k2 + (32.0 / 9.0) * k3),
        eos,
        params
    );
    const State k5 = rhs(
        r + h * (8.0 / 9.0),
        y + h * ((19372.0 / 6561.0) * k1 + (-25360.0 / 2187.0) * k2
            + (64448.0 / 6561.0) * k3 + (-212.0 / 729.0) * k4),
        eos,
        params
    );
    const State k6 = rhs(
        r + h,
        y + h * ((9017.0 / 3168.0) * k1 + (-355.0 / 33.0) * k2
            + (46732.0 / 5247.0) * k3 + (49.0 / 176.0) * k4 + (-5103.0 / 18656.0) * k5),
        eos,
        params
    );
    const State k7 = rhs(
        r + h,
        y + h * ((35.0 / 384.0) * k1 + (500.0 / 1113.0) * k3
            + (125.0 / 192.0) * k4 + (-2187.0 / 6784.0) * k5 + (11.0 / 84.0) * k6),
        eos,
        params
    );

    const State y5 = y + h * ((35.0 / 384.0) * k1 + (500.0 / 1113.0) * k3
        + (125.0 / 192.0) * k4 + (-2187.0 / 6784.0) * k5 + (11.0 / 84.0) * k6);
    const State y4 = y + h * ((5179.0 / 57600.0) * k1 + (7571.0 / 16695.0) * k3
        + (393.0 / 640.0) * k4 + (-92097.0 / 339200.0) * k5
        + (187.0 / 2100.0) * k6 + (1.0 / 40.0) * k7);

    return {y5, y5 - y4};
}

double normalized_error(const State& y, const State& ynext, const State& error,
                        double atol, double rtol, bool include_tidal) {
    const double scale_m = atol + rtol * std::max(std::abs(y.m), std::abs(ynext.m));
    const double scale_p = atol + rtol * std::max(std::abs(y.p), std::abs(ynext.p));
    const double em = error.m / scale_m;
    const double ep = error.p / scale_p;
    if (!include_tidal) {
        return std::sqrt(0.5 * (em * em + ep * ep));
    }
    const double scale_y = atol + rtol * std::max(std::abs(y.y), std::abs(ynext.y));
    const double ey = error.y / scale_y;
    return std::sqrt((em * em + ep * ep + ey * ey) / 3.0);
}

} // namespace

double weyl_suppression(double rho, const BraneParams& params) {
    if (params.closure != WeylClosure::V1Suppressed || params.rho0_geom <= 0.0) {
        return 1.0;
    }
    return 1.0 / (1.0 + std::pow(params.rho0_geom / rho, params.n_suppress));
}

double rho_weyl(double rho, const BraneParams& params) {
    return params.alpha_U * rho * weyl_suppression(rho, params);
}

double rho_eff_weyl(double rho, const BraneParams& params) {
    double rho_loc = rho;
    if (!std::isinf(params.lambda_brane)) {
        rho_loc = rho + rho * rho / (2.0 * params.lambda_brane);
    }
    return rho_loc + rho_weyl(rho, params);
}

double p_eff_weyl(double p, double rho, const BraneParams& params) {
    double p_loc = p;
    if (!std::isinf(params.lambda_brane)) {
        p_loc = p + p * rho / params.lambda_brane + rho * rho / (2.0 * params.lambda_brane);
    }
    // V2 ties the Weyl pressure to the matter pressure, so it vanishes with p at
    // the surface; the other closures tie it to the Weyl density.
    const double p_U = params.closure == WeylClosure::V2PressureTracking
        ? params.w_U * p
        : params.w_U * rho_weyl(rho, params);
    return p_loc + p_U;
}

double cs2_eff_weyl(double rho, double p, double cs2, const BraneParams& params) {
    if (params.closure == WeylClosure::V2PressureTracking) {
        // rho_tot = (1+a) rho + rho^2/2L ; p_tot = p(1+w) + p rho/L + rho^2/2L
        double drho_tot = 1.0 + params.alpha_U;
        double dp_tot = cs2 * (1.0 + params.w_U);
        if (!std::isinf(params.lambda_brane)) {
            drho_tot += rho / params.lambda_brane;
            dp_tot += (cs2 * rho + p) / params.lambda_brane + rho / params.lambda_brane;
        }
        return dp_tot / drho_tot;
    }
    if (params.closure == WeylClosure::V1Suppressed && params.rho0_geom > 0.0) {
        // rho_U = a rho f(rho); d(rho_U)/drho = a (f + rho f')
        const double x = params.rho0_geom / rho;
        const double xn = std::pow(x, params.n_suppress);
        const double f = 1.0 / (1.0 + xn);
        const double dfdrho = params.n_suppress * xn / (rho * (1.0 + xn) * (1.0 + xn));
        const double drhoU = params.alpha_U * (f + rho * dfdrho);
        double drho_tot = 1.0 + drhoU;
        double dp_tot = cs2 + params.w_U * drhoU;
        if (!std::isinf(params.lambda_brane)) {
            drho_tot += rho / params.lambda_brane;
            dp_tot += (cs2 * rho + p) / params.lambda_brane + rho / params.lambda_brane;
        }
        return dp_tot / drho_tot;
    }
    //   rho_tot = (1 + alpha) rho + rho^2/(2 lambda)
    //   p_r,tot = p + p rho / lambda + rho^2/(2 lambda) + w alpha rho
    // Both depend on r only through rho, so
    //   d p_r,tot / d rho_tot = (d p_r,tot/d rho) / (d rho_tot/d rho).
    double drho_tot = 1.0 + params.alpha_U;
    double dp_tot = cs2 + params.w_U * params.alpha_U;
    if (!std::isinf(params.lambda_brane)) {
        drho_tot += rho / params.lambda_brane;
        dp_tot += (cs2 * rho + p) / params.lambda_brane + rho / params.lambda_brane;
    }
    return dp_tot / drho_tot;
}

Star integrate_star(double rho_c, const EOS& eos, const BraneParams& params, double r_max, double dr, double r0, std::vector<ProfilePoint>* profile) {
    const double p_c = eos.p_of_rho(rho_c);
    if (!std::isfinite(p_c) || p_c <= 0.0) {
        throw std::runtime_error("Pressao central invalida.");
    }

    const IntegratorConfig& icfg = integrator_config();
    const double rho_tot_c = rho_eff_weyl(rho_c, params);
    State y{4.0 * pi * r0 * r0 * r0 * rho_tot_c / 3.0, p_c};
    double r = r0;
    // With an analytic atmosphere p_min() is 0, and the event limiter below makes p
    // decay geometrically without ever reaching it, so the integration would run to
    // r_max. The surface therefore needs a declared floor relative to p_c. For a
    // polytrope p ~ (r_s - r)^(Gamma/(Gamma-1)), so the radius left beyond the floor
    // scales as p_stop^((Gamma-1)/Gamma) and a small floor costs few steps.
    const double p_stop_frac = icfg.p_stop_frac;
    const double p_surface = eos.p_min() > 0.0 ? eos.p_min() * 1.0001
                                               : p_stop_frac * p_c;
    const double rtol = icfg.rtol;
    const double sigma_surface = icfg.sigma_surface;
    const double atol = icfg.atol;
    constexpr double h_min = 1.0e-6;

    auto record = [&](double rr, const State& yy) {
        if (!profile) {
            return;
        }
        const double rho = eos.rho_of_p(yy.p);
        if (!std::isfinite(rho)) {
            return;
        }
        const double cs2 = eos.cs2_of_rho(rho);
        ProfilePoint pt;
        pt.r = rr;
        pt.m = yy.m;
        pt.p = yy.p;
        pt.rho = rho;
        pt.rho_tot = rho_eff_weyl(rho, params);
        pt.p_r_tot = p_eff_weyl(yy.p, rho, params);
        pt.cs2 = cs2;
        pt.cs2_eff = cs2_eff_weyl(rho, yy.p, cs2, params);
        profile->push_back(pt);
    };

    double prev_r = r;
    State prev_y = y;
    double h = std::min(dr, 1.0e-3);
    record(r, y);
    int n_internal_steps = 0;
    constexpr int max_internal_steps = 200000;
    while (r < r_max && y.p > p_surface && std::isfinite(y.m) && std::isfinite(y.p)) {
        if (++n_internal_steps > max_internal_steps) {
            throw std::runtime_error("Integracao TOV excedeu o limite de passos internos.");
        }
        h = std::min({h, dr, r_max - r});

        // Event-based step limiter. Near the surface p falls below atol, so the
        // relative error test compares against a floor larger than the quantity
        // itself: the adaptive control silently switches off exactly where the
        // radius is decided, steps stay pinned at the dr cap, and the integration
        // overshoots into p < 0 before any bisection can run (measured: p =
        // -2.7e-13 at the last layer). Capping h so that a single step cannot
        // consume more than a fraction sigma of the remaining pressure keeps the
        // surface crossing bracketed. Far from the surface p/|dp/dr| is large and
        // this does not bind. See docs/PHASE_C1.md.
        if (y.p > 0.0 && icfg.use_surface_limiter) {
            const State d = rhs(r, y, eos, params);
            if (d.p < 0.0) {
                h = std::min(h, sigma_surface * y.p / (-d.p));
            }
        }
        const RKStep trial = rk45_step(r, y, h, eos, params);
        const double err = normalized_error(y, trial.y5, trial.error, atol, rtol,
                                            icfg.tidal_in_error_control);

        if (err <= 1.0 || h <= h_min) {
            prev_r = r;
            prev_y = y;
            y = trial.y5;
            r += h;

            const double factor = err > 0.0
                ? std::clamp(0.9 * std::pow(1.0 / err, 0.2), 0.2, 5.0)
                : 5.0;
            h = std::min(dr, h * factor);
        } else {
            const double factor = std::clamp(0.9 * std::pow(1.0 / err, 0.25), 0.1, 0.5);
            h *= factor;
            if (h < h_min) {
                h = h_min;
            }
            continue;
        }

        record(r, y);

        if (r <= 2.0 * y.m || !std::isfinite(y.m) || !std::isfinite(y.p)) {
            break;
        }
    }

    double R = r;
    double M = y.m;
    double Y_R = y.y;
    if (y.p <= p_surface && prev_y.p > p_surface) {
        double left_r = prev_r;
        State left_y = prev_y;
        double right_r = r;
        State right_y = y;

        for (int i = 0; i < 32; ++i) {
            const double mid_r = 0.5 * (left_r + right_r);
            const double mid_h = mid_r - left_r;
            const State mid_y = rk45_step(left_r, left_y, mid_h, eos, params).y5;
            if (mid_y.p > p_surface) {
                left_r = mid_r;
                left_y = mid_y;
            } else {
                right_r = mid_r;
                right_y = mid_y;
            }
        }

        R = right_r;
        M = right_y.m;
        Y_R = right_y.y;
    }

    Star star;
    star.y_R = Y_R;
    star.R = R;
    star.M_km = M;
    star.M_sun = M / M_sun_km;
    star.compactness = R > 0.0 ? M / R : std::numeric_limits<double>::quiet_NaN();
    star.rho_c = rho_c;
    star.rho_c_cgs = rho_geom_km2_to_cgs(rho_c);
    star.p_c = p_c;
    star.params = params;
    return star;
}

Sequence mass_radius_sequence(
    const std::vector<double>& rho_c_values,
    const EOS& eos,
    const BraneParams& params,
    const std::string& label,
    double r_max,
    double dr
) {
    Sequence seq;
    seq.label = label;
    for (double rho_c : rho_c_values) {
        try {
            Star star = integrate_star(rho_c, eos, params, r_max, dr);
            if (std::isfinite(star.R) && std::isfinite(star.M_sun) && star.R > 0.0 && star.M_sun > 0.0) {
                seq.stars.push_back(star);
            }
        } catch (const std::exception&) {
        }
    }
    return seq;
}

Diagnostics sequence_diagnostics(const Sequence& seq) {
    Diagnostics d;
    if (seq.stars.empty()) {
        return d;
    }

    auto imax_it = std::max_element(seq.stars.begin(), seq.stars.end(), [](const Star& a, const Star& b) {
        return a.M_sun < b.M_sun;
    });
    const auto imax = static_cast<std::size_t>(std::distance(seq.stars.begin(), imax_it));
    d.valid = true;
    d.Mmax = seq.stars[imax].M_sun;
    d.Rmax = seq.stars[imax].R;
    d.rhoc_Mmax = seq.stars[imax].rho_c_cgs;
    d.Cmax = seq.stars[imax].compactness;

    std::vector<std::pair<double, double>> mr;
    for (std::size_t i = 0; i <= imax; ++i) {
        mr.emplace_back(seq.stars[i].M_sun, seq.stars[i].R);
    }
    std::sort(mr.begin(), mr.end());

    constexpr double target = 1.4;
    if (mr.size() >= 2 && target >= mr.front().first && target <= mr.back().first) {
        for (std::size_t i = 0; i + 1 < mr.size(); ++i) {
            if (target >= mr[i].first && target <= mr[i + 1].first) {
                const double denom = mr[i + 1].first - mr[i].first;
                const double t = denom != 0.0 ? (target - mr[i].first) / denom : 0.0;
                d.R14 = mr[i].second + t * (mr[i + 1].second - mr[i].second);
                break;
            }
        }
    }

    return d;
}

std::vector<double> logspace(double log10_min, double log10_max, int n) {
    std::vector<double> values;
    values.reserve(static_cast<std::size_t>(n));
    if (n <= 1) {
        values.push_back(std::pow(10.0, log10_min));
        return values;
    }
    for (int i = 0; i < n; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(n - 1);
        values.push_back(std::pow(10.0, log10_min + t * (log10_max - log10_min)));
    }
    return values;
}

void write_sequence_dat(const Sequence& seq, const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Nao foi possivel escrever: " + path);
    }
    out << "# label " << seq.label << "\n";
    out << "# rho_c_cgs rho_c_geom p_c_geom R_km M_km M_sun compactness lambda alpha_U w_U y_R\n";
    out << std::setprecision(16);
    for (const auto& s : seq.stars) {
        out << s.rho_c_cgs << ' ' << s.rho_c << ' ' << s.p_c << ' ' << s.R << ' ' << s.M_km << ' '
            << s.M_sun << ' ' << s.compactness << ' ' << s.params.lambda_brane << ' ' << s.params.alpha_U
            << ' ' << s.params.w_U << ' ' << s.y_R << '\n';
    }
}

void write_stable_sequence_dat(const Sequence& seq, const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Nao foi possivel escrever: " + path);
    }
    out << "# label " << seq.label << " stable_branch_until_Mmax\n";
    out << "# rho_c_cgs rho_c_geom p_c_geom R_km M_km M_sun compactness lambda alpha_U w_U y_R\n";
    out << std::setprecision(16);

    if (seq.stars.empty()) {
        return;
    }

    auto imax_it = std::max_element(seq.stars.begin(), seq.stars.end(), [](const Star& a, const Star& b) {
        return a.M_sun < b.M_sun;
    });
    const auto imax = static_cast<std::size_t>(std::distance(seq.stars.begin(), imax_it));

    for (std::size_t i = 0; i <= imax; ++i) {
        const auto& s = seq.stars[i];
        out << s.rho_c_cgs << ' ' << s.rho_c << ' ' << s.p_c << ' ' << s.R << ' ' << s.M_km << ' '
            << s.M_sun << ' ' << s.compactness << ' ' << s.params.lambda_brane << ' ' << s.params.alpha_U
            << ' ' << s.params.w_U << ' ' << s.y_R << '\n';
    }
}

} // namespace brana
