#include "eos.hpp"

#include "constants.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace brana {

MITBagEOS::MITBagEOS(double B_geom, double a) : B_(B_geom), a_(a) {}

double MITBagEOS::p_of_rho(double rho) const {
    return a_ * (rho - 4.0 * B_);
}

double MITBagEOS::rho_of_p(double p) const {
    return p / a_ + 4.0 * B_;
}

double MITBagEOS::cs2_of_rho(double) const { return a_; }

double MITBagEOS::rho_min() const { return rho_surface(); }
double MITBagEOS::rho_max() const { return std::numeric_limits<double>::infinity(); }
double MITBagEOS::p_min() const { return 0.0; }
double MITBagEOS::p_max() const { return std::numeric_limits<double>::infinity(); }
double MITBagEOS::rho_surface() const { return 4.0 * B_; }

std::string MITBagEOS::name() const {
    std::ostringstream os;
    os << "MIT bag: a=" << a_ << ", B=" << B_ << " km^-2";
    return os.str();
}

TabulatedEOS::TabulatedEOS(std::vector<double> rho_geom, std::vector<double> p_geom, std::string name)
    : name_(std::move(name)) {
    if (rho_geom.size() != p_geom.size() || rho_geom.size() < 2) {
        throw std::runtime_error("EOS tabulada precisa de colunas rho e P com pelo menos dois pontos.");
    }

    std::vector<std::pair<double, double>> rows;
    rows.reserve(rho_geom.size());
    for (std::size_t i = 0; i < rho_geom.size(); ++i) {
        if (std::isfinite(rho_geom[i]) && std::isfinite(p_geom[i]) && rho_geom[i] > 0.0 && p_geom[i] > 0.0) {
            rows.emplace_back(rho_geom[i], p_geom[i]);
        }
    }

    std::sort(rows.begin(), rows.end(), [](const auto& a, const auto& b) { return a.first < b.first; });

    double last_rho = -1.0;
    double last_p = -1.0;
    for (const auto& row : rows) {
        if (row.first <= last_rho) {
            continue;
        }
        if (!p_.empty() && row.second <= last_p) {
            continue;
        }
        rho_.push_back(row.first);
        p_.push_back(row.second);
        last_rho = row.first;
        last_p = row.second;
    }

    if (rho_.size() < 2) {
        throw std::runtime_error("EOS tabulada ficou com menos de dois pontos monotônicos.");
    }

    logrho_.reserve(rho_.size());
    logp_.reserve(p_.size());
    for (std::size_t i = 0; i < rho_.size(); ++i) {
        logrho_.push_back(std::log(rho_[i]));
        logp_.push_back(std::log(p_[i]));
    }

    logp_of_logrho_ = Pchip(logrho_, logp_);
    logrho_of_logp_ = Pchip(logp_, logrho_);
}

// Automatic rule for rho_join: the lowest tabulated density whose local logarithmic
// slope stays inside a well-behaved band over [rho, 10 rho]. Six of the seven EoSs
// land on their own table minimum (nothing is replaced); only SLy, which tabulates
// down to 1 g/cm^3 with Gamma reaching ~5, has its tail actually discarded.
void TabulatedEOS::enable_atmosphere(double rho_join_geom) {
    constexpr double GAMMA_LO = 1.2;
    constexpr double GAMMA_HI = 2.0;

    if (rho_join_geom <= 0.0) {
        rho_join_geom = rho_.front();
        for (std::size_t i = 0; i < rho_.size(); ++i) {
            const double lo = rho_[i];
            const double hi = 10.0 * lo;
            if (hi > rho_.back()) {
                break;
            }
            double gmin = std::numeric_limits<double>::infinity();
            double gmax = -std::numeric_limits<double>::infinity();
            for (std::size_t j = i; j < rho_.size() && rho_[j] <= hi; ++j) {
                const double g = logp_of_logrho_.derivative(logrho_[j]);
                gmin = std::min(gmin, g);
                gmax = std::max(gmax, g);
            }
            if (gmin >= GAMMA_LO && gmax <= GAMMA_HI) {
                rho_join_geom = lo;
                break;
            }
        }
    }
    rho_join_geom = std::clamp(rho_join_geom, rho_.front(), rho_.back());

    rho_join_ = rho_join_geom;
    const double logrj = std::log(rho_join_);
    gamma_atm_ = logp_of_logrho_.derivative(logrj);
    if (!(gamma_atm_ > 1.0)) {
        throw std::runtime_error("Atmosfera exige Gamma > 1 em rho_join; obtido "
                                 + std::to_string(gamma_atm_));
    }
    p_join_ = std::exp(logp_of_logrho_(logrj));
    k_atm_ = p_join_ / std::pow(rho_join_, gamma_atm_);
    atm_ = true;
}

TabulatedEOS TabulatedEOS::from_file_cgs(const std::string& path, double rho_join_cgs) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("Nao foi possivel abrir EOS: " + path);
    }

    std::vector<double> rho;
    std::vector<double> pressure;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::istringstream is(line);
        double rho_cgs = 0.0;
        double p_cgs = 0.0;
        if (is >> rho_cgs >> p_cgs) {
            rho.push_back(rho_mass_cgs_to_geom_km2(rho_cgs));
            pressure.push_back(pressure_cgs_to_geom_km2(p_cgs));
        }
    }

    TabulatedEOS eos(std::move(rho), std::move(pressure), "EoS tabulada: " + path);
    if (rho_join_cgs >= 0.0) {
        eos.enable_atmosphere(rho_join_cgs > 0.0 ? rho_mass_cgs_to_geom_km2(rho_join_cgs) : 0.0);
    }
    return eos;
}

double TabulatedEOS::p_of_rho(double rho) const {
    if (atm_ && rho < rho_join_) {
        return rho > 0.0 ? k_atm_ * std::pow(rho, gamma_atm_) : 0.0;
    }
    return std::exp(logp_of_logrho_(std::log(rho)));
}

double TabulatedEOS::rho_of_p(double p) const {
    if (atm_ && p < p_join_) {
        return p > 0.0 ? std::pow(p / k_atm_, 1.0 / gamma_atm_) : 0.0;
    }
    return std::exp(logrho_of_logp_(std::log(p)));
}

// With log p interpolated against log rho, dp/drho = (p/rho) * dlogp/dlogrho.
double TabulatedEOS::cs2_of_rho(double rho) const {
    if (atm_ && rho < rho_join_) {
        // dp/drho = Gamma K rho^(Gamma-1) = Gamma p / rho, continuous at the join
        // because Gamma is the local slope there.
        return rho > 0.0 ? gamma_atm_ * k_atm_ * std::pow(rho, gamma_atm_ - 1.0) : 0.0;
    }
    const double logrho = std::log(rho);
    const double slope = logp_of_logrho_.derivative(logrho);
    return std::exp(logp_of_logrho_(logrho)) / rho * slope;
}

double TabulatedEOS::rho_min() const { return atm_ ? 0.0 : rho_.front(); }
double TabulatedEOS::rho_max() const { return rho_.back(); }
double TabulatedEOS::p_min() const { return atm_ ? 0.0 : p_.front(); }
double TabulatedEOS::p_max() const { return p_.back(); }
std::string TabulatedEOS::name() const { return name_; }

} // namespace brana
