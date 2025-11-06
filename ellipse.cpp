#include "ellipse.h"

#include <cmath>
#include <stdexcept>
#include <algorithm>

Ellipse::Ellipse(double xc_i, double yc_i, double a_i, double b_i, double theta_i){
        double final_a = a_i;
        double final_b = b_i;
        double final_theta = theta_i;

        // enforce a >= b
        if (final_a < final_b) {
            std::swap(final_a, final_b);
            final_theta += pi / 2.0;
        }

        // enforce theta in [0,2*pi)
        final_theta = std::fmod(final_theta, two_pi);
        if (final_theta < 0.0) {
            final_theta += two_pi;
        }

        //check if ellipse is inside the square with margin h
        check_bounds(xc_i, yc_i, final_a, final_b, final_theta);

        xc = xc_i;
        yc = yc_i;
        a = final_a;
        b = final_b;
        theta = final_theta;
}

void Ellipse::meshify(){
        WorkingManifold RR2WM;
        Manifold &RR2 = RR2WM.current;
        Function xy = RR2.coordinates();
        Function x = xy[0], y = xy[1];

        double sin_t = std::sin(theta), cos_t = std::cos(theta);
        if (! representation.has_value()){
            Function left_numerator  = (x - xc) * cos_t + (y - yc) * sin_t;
            Function right_numerator = (x - xc) * sin_t - (y - yc) * cos_t;
            representation.emplace(RR2.implicit((left_numerator * left_numerator)/(a*a) + (right_numerator * right_numerator)/(b*b) == 1));
        }

        
        representation->set_as_working_manifold();
        Cell A(tag::vertex, tag::of_coords, {xc + a * cos_t , yc + a * sin_t}); //extremidade equivalente a (a,0)
        std::vector<double> direction = {sin_t, - cos_t}; //direção equivalente a (0,-1)
        mesh.emplace(Mesh::Build(tag::frontal).entire_manifold().start_at(A).towards(direction).desired_length(h));
}

void Ellipse::check_bounds(double xc_i, double yc_i, double a_i, double b_i, double theta_i) const {
    double sin_t = std::sin(theta_i), cos_t = std::cos(theta_i);

    double min_height = h + std::sqrt(a_i * a_i * cos_t * cos_t + b_i * b_i * sin_t * sin_t);
    double min_width = h + std::sqrt(a_i * a_i * sin_t * sin_t + b_i * b_i * cos_t * cos_t);

    if ((yc_i < min_height) || (xc_i < min_width) || (yc_i > 1 - min_height) || (xc_i > 1 - min_width) ){
        throw std::invalid_argument("Ellipse does not fit in the unit square.");
    }
}