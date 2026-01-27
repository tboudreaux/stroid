#pragma once

#include <functional>
#include <variant>

namespace stroid::core {
    using SpacingFunction = std::function<double(double)>

    struct LinearSpacing {
        double operator()(double xi) const {return xi;}
    };

    struct LogarithmicSpacing {
        double base = 10.0;
        double operator()(double xi) const {
            return (std::pow(base, xi) - 1) / (base - 1.0);
        }
    };

    struct GeometricSpacing {
        double ratio = 1.2;
        double operator()(double xi) const {
            return (std::pow(ratio, xi) - 1) / (ratio - 1.0);
        }
    };

    using SpacingStrategy = std::variant<LinearSpacing, LogarithmicSpacing, GeometricSpacing, SpacingFunction>;
}