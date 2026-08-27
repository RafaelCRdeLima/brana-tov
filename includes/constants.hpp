#pragma once

#include <cmath>

namespace brana {

constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double G_cgs = 6.67430e-8;
constexpr double c_cgs = 2.99792458e10;
constexpr double M_sun_km = 1.47662504;

inline double rho_mass_cgs_to_geom_km2(double rho_cgs) {
    return rho_cgs * G_cgs / (c_cgs * c_cgs) * 1.0e10;
}

inline double pressure_cgs_to_geom_km2(double p_cgs) {
    return p_cgs * G_cgs / std::pow(c_cgs, 4) * 1.0e10;
}

inline double rho_geom_km2_to_cgs(double rho_geom) {
    return rho_geom * c_cgs * c_cgs / G_cgs / 1.0e10;
}

inline double pressure_geom_km2_to_cgs(double p_geom) {
    return p_geom * std::pow(c_cgs, 4) / G_cgs / 1.0e10;
}

inline double MeVfm3_to_dyn_cm2(double x) {
    return x * 1.602176634e33;
}

} // namespace brana
