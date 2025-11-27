#ifndef DOMAIN_CONFIG_H
#define DOMAIN_CONFIG_H

#include "maniFEM.h"
#include "constants.h"

class DomainConfig {
    public:
        const double x_max, y_max;
        const double MW_x, ME_x;
        const double h;
        const size_t num_ellipses;

        DomainConfig() : x_max(1), y_max(1), MW_x(0.3), ME_x(0.7), h(0.02), num_ellipses(10) {};

        DomainConfig(const double x_m,const double y_m, const double MW, const double ME, const double h_size, const size_t n_ellipses) :
        x_max(x_m), y_max(y_m), MW_x(MW), ME_x(ME), h(h_size), num_ellipses(n_ellipses) {
            if (x_max < 0 || y_max < 0 || //retângulo inválido
            MW_x > ME_x || ME_x > x_max || // MW à direita de ME ou ME à direita de y_max
            h < 0){
                throw std::invalid_argument("Invalid Domain Configuration");
            }
        }

        int n_segments(double l){ return std::ceil(l/h); }
};

#endif