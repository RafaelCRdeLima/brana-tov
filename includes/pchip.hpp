#pragma once

#include <stdexcept>
#include <vector>

namespace brana {

class Pchip {
public:
    Pchip() = default;
    Pchip(std::vector<double> x, std::vector<double> y);

    double operator()(double xq) const;
    double derivative(double xq) const;
    bool in_range(double xq) const;

private:
    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> d_;
};

} // namespace brana
