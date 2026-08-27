#include "pchip.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace brana {
namespace {

double endpoint_slope(double h0, double h1, double delta0, double delta1) {
    double d = ((2.0 * h0 + h1) * delta0 - h0 * delta1) / (h0 + h1);
    if (d * delta0 <= 0.0) {
        return 0.0;
    }
    if (delta0 * delta1 < 0.0 && std::abs(d) > std::abs(3.0 * delta0)) {
        return 3.0 * delta0;
    }
    return d;
}

} // namespace

Pchip::Pchip(std::vector<double> x, std::vector<double> y) : x_(std::move(x)), y_(std::move(y)) {
    const std::size_t n = x_.size();
    if (n != y_.size() || n < 2) {
        throw std::runtime_error("PCHIP precisa de pelo menos dois pontos.");
    }

    std::vector<double> h(n - 1);
    std::vector<double> delta(n - 1);

    for (std::size_t i = 0; i + 1 < n; ++i) {
        h[i] = x_[i + 1] - x_[i];
        if (h[i] <= 0.0) {
            throw std::runtime_error("PCHIP exige eixo x estritamente crescente.");
        }
        delta[i] = (y_[i + 1] - y_[i]) / h[i];
    }

    d_.assign(n, 0.0);
    if (n == 2) {
        d_[0] = delta[0];
        d_[1] = delta[0];
        return;
    }

    d_[0] = endpoint_slope(h[0], h[1], delta[0], delta[1]);
    d_[n - 1] = endpoint_slope(h[n - 2], h[n - 3], delta[n - 2], delta[n - 3]);

    for (std::size_t k = 1; k + 1 < n; ++k) {
        if (delta[k - 1] == 0.0 || delta[k] == 0.0 || delta[k - 1] * delta[k] < 0.0) {
            d_[k] = 0.0;
        } else {
            const double w1 = 2.0 * h[k] + h[k - 1];
            const double w2 = h[k] + 2.0 * h[k - 1];
            d_[k] = (w1 + w2) / (w1 / delta[k - 1] + w2 / delta[k]);
        }
    }
}

bool Pchip::in_range(double xq) const {
    return !x_.empty() && xq >= x_.front() && xq <= x_.back();
}

double Pchip::operator()(double xq) const {
    if (!in_range(xq)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    auto upper = std::upper_bound(x_.begin(), x_.end(), xq);
    std::size_t i = 0;
    if (upper == x_.begin()) {
        i = 0;
    } else if (upper == x_.end()) {
        i = x_.size() - 2;
    } else {
        i = static_cast<std::size_t>(std::distance(x_.begin(), upper) - 1);
    }

    const double h = x_[i + 1] - x_[i];
    const double t = (xq - x_[i]) / h;
    const double t2 = t * t;
    const double t3 = t2 * t;

    const double h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    const double h10 = t3 - 2.0 * t2 + t;
    const double h01 = -2.0 * t3 + 3.0 * t2;
    const double h11 = t3 - t2;

    return h00 * y_[i] + h10 * h * d_[i] + h01 * y_[i + 1] + h11 * h * d_[i + 1];
}

double Pchip::derivative(double xq) const {
    if (!in_range(xq)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    auto upper = std::upper_bound(x_.begin(), x_.end(), xq);
    std::size_t i = 0;
    if (upper == x_.begin()) {
        i = 0;
    } else if (upper == x_.end()) {
        i = x_.size() - 2;
    } else {
        i = static_cast<std::size_t>(std::distance(x_.begin(), upper) - 1);
    }

    const double h = x_[i + 1] - x_[i];
    const double t = (xq - x_[i]) / h;
    const double t2 = t * t;

    // d/dt of the cubic Hermite basis, divided by h.
    const double dh00 = 6.0 * t2 - 6.0 * t;
    const double dh10 = 3.0 * t2 - 4.0 * t + 1.0;
    const double dh01 = -6.0 * t2 + 6.0 * t;
    const double dh11 = 3.0 * t2 - 2.0 * t;

    return (dh00 * y_[i] + dh01 * y_[i + 1]) / h + dh10 * d_[i] + dh11 * d_[i + 1];
}

} // namespace brana
