#ifndef ELLIPSE_H
#define ELLIPSE_H

#include "maniFEM.h"
#include "constants.h"
#include <optional>
#include <numbers>

// A(x-xc)^2 + 2B(x-xc)(y-yc) + C(y-yc)^2 = 1
struct QuadraticForm {
    double A, B, C;
};

class WorkingManifold{
    public:
        WorkingManifold() = default;
        ~WorkingManifold() {current.set_as_working_manifold();}

        WorkingManifold(const WorkingManifold&) = delete;
        WorkingManifold& operator=(const WorkingManifold&) = delete;

        Manifold &current = Manifold::working;
};

class Ellipse {
    public:
        double xc, yc, a, b, theta, height, width; 
        std::optional<Manifold> representation;
        std::optional<Mesh> mesh;
        std::optional<QuadraticForm> q_form;

        Ellipse(double xc_i, double yc_i, double a_i, double b_i, double theta_i);

        void meshify();

        Mesh get_mesh() {
            if (! mesh.has_value()){
                meshify();
            }
            return mesh.value();
        }

        double area(){return pi*a*b;}
        //pela rota da forma quadrática pode ser pi*det(M)

        QuadraticForm& quadratic_form(){
            if (! q_form.has_value()){
                q_form = calculate_quadratic_form();
            }
            return q_form.value();
        }
        
    private:
        void check_bounds(double xc_i, double yc_i, double height_i, double width_i) const;
        QuadraticForm calculate_quadratic_form();
};

class EllipseBundle{
    public:
        std::vector<Ellipse> bundle;
        
        EllipseBundle(){bundle.reserve(num_ellipses);}

        void add_ellipse(const Ellipse &new_ellipse){
            check_intersections(new_ellipse);
            bundle.push_back(new_ellipse);
        }

    private:
        void check_intersections(const Ellipse &new_ellipse){
            for (auto & ellipse : bundle){
                if (intersects(new_ellipse, ellipse)){ 
                    throw std::invalid_argument("Ellipses might intersect");
                };
            }
        };

        bool intersects(const Ellipse &Ellipse1, const Ellipse &Ellipse2);
        bool robust_intersect(const Ellipse &Ellipse1, const Ellipse &Ellipse2); //not yet implemented
};
#endif