#pragma once

#include "pchip.hpp"

#include <memory>
#include <string>
#include <vector>

namespace brana {

class EOS {
public:
    virtual ~EOS() = default;
    virtual double p_of_rho(double rho) const = 0;
    virtual double rho_of_p(double p) const = 0;
    // Squared sound speed of the *matter* sector, c_s^2 = dp/drho (units c=1).
    virtual double cs2_of_rho(double rho) const = 0;
    virtual double rho_min() const = 0;
    virtual double rho_max() const = 0;
    virtual double p_min() const = 0;
    virtual double p_max() const = 0;
    virtual std::string name() const = 0;
};

class MITBagEOS final : public EOS {
public:
    MITBagEOS(double B_geom, double a = 1.0 / 3.0);

    double p_of_rho(double rho) const override;
    double rho_of_p(double p) const override;
    double cs2_of_rho(double rho) const override;
    double rho_min() const override;
    double rho_max() const override;
    double p_min() const override;
    double p_max() const override;
    std::string name() const override;

    double rho_surface() const;

private:
    double B_;
    double a_;
};

class TabulatedEOS final : public EOS {
public:
    TabulatedEOS(std::vector<double> rho_geom, std::vector<double> p_geom, std::string name);
    // rho_join_cgs < 0 : no atmosphere (tabulated floor is the surface, legacy behaviour)
    //              = 0 : automatic rule, see enable_atmosphere()
    //              > 0 : explicit join density, for the sensitivity sweep
    static TabulatedEOS from_file_cgs(const std::string& path, double rho_join_cgs = -1.0);

    // Replaces the tabulated tail below rho_join by an analytic polytrope
    // p = K rho^Gamma integrable to p = 0. Gamma is the LOCAL logarithmic slope at
    // rho_join and K is anchored on p(rho_join), so p and dp/drho are both continuous
    // there by construction -- a least-squares Gamma over a band would match p but
    // not dp/drho. Choosing rho_join inside the well-behaved region is what keeps
    // Gamma physical: SLy's tabulated tail reaches Gamma ~ 5, which is an artefact of
    // its low-density fit and must be replaced, not continued. See docs/PHASE_C1.md.
    void enable_atmosphere(double rho_join_geom = 0.0);
    bool has_atmosphere() const { return atm_; }
    double atmosphere_rho_join() const { return rho_join_; }
    double atmosphere_gamma() const { return gamma_atm_; }

    double p_of_rho(double rho) const override;
    double rho_of_p(double p) const override;
    double cs2_of_rho(double rho) const override;
    double rho_min() const override;
    double rho_max() const override;
    double p_min() const override;
    double p_max() const override;
    std::string name() const override;

private:
    std::vector<double> rho_;
    std::vector<double> p_;
    std::vector<double> logrho_;
    std::vector<double> logp_;
    Pchip logp_of_logrho_;
    Pchip logrho_of_logp_;
    std::string name_;
    bool atm_ = false;
    double rho_join_ = 0.0;
    double p_join_ = 0.0;
    double gamma_atm_ = 0.0;
    double k_atm_ = 0.0;
};

} // namespace brana
