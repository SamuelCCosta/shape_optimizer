#include "ellipse_old.h"

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

    // enforce theta in [0,pi)
    final_theta = std::fmod(final_theta, pi);
    if (final_theta < 0.0) {
        final_theta += pi;
    }

    double sin_t = std::sin(theta_i), cos_t = std::cos(theta_i);
    double width_i = std::sqrt(a_i * a_i * cos_t * cos_t + b_i * b_i * sin_t * sin_t);
    // este pode ser sqrt(A/det(M)) se for pela rota da forma quadrática
    double height_i = std::sqrt(a_i * a_i * sin_t * sin_t + b_i * b_i * cos_t * cos_t);
    // este pode ser sqrt(C/det(M)) se for pela rota da forma quadrática
    // Fonte para os valores em formas quadráticas: https://www.geometrictools.com/Documentation/RobustIntersectionOfEllipses.pdf

    //check if ellipse is inside the square with margin h
    check_bounds(xc_i, yc_i, width_i, height_i);

    xc = xc_i;
    yc = yc_i;
    a = final_a;
    b = final_b;
    theta = final_theta;
    width = width_i;
    height = height_i;
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

void Ellipse::check_bounds(double xc_i, double yc_i, double height_i, double width_i) const {
    double horizontal_margin = h + width_i;
    double vertical_margin = h + height_i;

    if ((yc_i < vertical_margin) || (xc_i < horizontal_margin) || (yc_i > 1 - vertical_margin) || (xc_i > 1 - horizontal_margin) ){
        throw std::invalid_argument("Ellipse does not fit in the unit square.");
    }
}

QuadraticForm Ellipse::calculate_quadratic_form(){
    double cos_t = std::cos(theta);
    double sin_t = std::sin(theta);
    double a2 = a * a;
    double b2 = b * b;

    double A = cos_t * cos_t / a2 + sin_t * sin_t / b2;
    //Ax^2 + 2Bxy + Cy^2
    double B = sin_t * cos_t * (1.0 / a2 - 1.0 / b2);
    double C = sin_t * sin_t / a2 + cos_t * cos_t / b2;

    return {A, B, C};    
}

bool EllipseBundle::intersects(const Ellipse &e1, const Ellipse &e2){
    //b = bottom, t = top, l = left, r = right
    double e1_l_x = e1.xc - e1.width; 
    double e1_b_y = e1.yc - e1.height;
    double e1_r_x = e1.xc + e1.width;
    double e1_t_y = e1.yc + e1.height;

    double e2_l_x = e2.xc - e2.width; 
    double e2_b_y = e2.yc - e2.height;
    double e2_r_x = e2.xc + e2.width;
    double e2_t_y = e2.yc + e2.height;

    if (e1_r_x < e2_l_x || //e1 à esquerda de e2
        e2_r_x < e1_l_x || //e1 à direita de e2
        e1_t_y < e2_b_y || //e1 abaixo de e2
        e2_t_y < e1_b_y)   //e1 acima de e2
    {return false;}
    
    return robust_intersect(e1,e2);
}

bool robust_intersect(const Ellipse &e1, const Ellipse &e2){
    return true;
}