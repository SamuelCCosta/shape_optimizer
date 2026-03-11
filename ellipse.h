#ifndef ELLIPSE_H
#define ELLIPSE_H

#include "maniFEM.h"
#include <map>
#include <string>
#include <numbers>
#include <Eigen/Dense>
#include <numbers>

constexpr double pi = std::numbers::pi;

// Class of ellipses with eccentricity < 0.93
class Ellipse {
    public:
        //(x-center)^T*q_form*(x-center) = 1
        // also known as (x-center, q_form*(x-center))
        
        Eigen::Vector2d center;
        Eigen::Matrix2d quadratic_form;
        Eigen::Matrix2d parametrize_matrix;
        Eigen::Vector2d bounds; //bounding half-heights
        
        // Enforces the ellipse invariants and eccentricity < 0.93
        Ellipse(double x, double y, double A, double B, double C);
        
        Mesh get_mesh(const double h) const;

        double area() const {return pi * 1 / std::sqrt(quadratic_form.determinant());}
        
        //https://www.geometrictools.com/Documentation/RobustIntersectionOfEllipses.pdf
        //bounding boxes half heights
        double width() const {return std::sqrt(quadratic_form(0,0)/quadratic_form.determinant());}

        double height() const {return std::sqrt(quadratic_form(1,1)/quadratic_form.determinant());}

        Eigen::Vector2d bounding_half() const {return Eigen::Vector2d(width(), height());}
        
        double evaluate_at(Eigen::Vector2d point) const {
            Eigen::Vector2d d = point - center;
            return d.dot(quadratic_form * d);
        }

        Eigen::Vector2d point_at(const double t) const { 
            Eigen::Vector2d point(std::cos(t), std::sin(t));
            return center + parametrize_matrix * point;
        }

        Eigen::Vector2d derivative_at(const double t) const {
            Eigen::Vector2d point(- std::sin(t), std::cos(t));
            return parametrize_matrix * point;
        }

        Eigen::Vector2d second_derivative_at(const double t) const {
            Eigen::Vector2d point(- std::cos(t),- std::sin(t));
            return parametrize_matrix * point;
        }

        // must be unitary vector
        Eigen::Vector2d point_at(Eigen::Vector2d point) const { 
            return center + parametrize_matrix * point;
        }

        // must be unitary vector
        Eigen::Vector2d derivative_at(Eigen::Vector2d point) const {
            return parametrize_matrix * point;
        }

        // must be unitary vector
        Eigen::Vector2d second_derivative_at(Eigen::Vector2d point) const {
            return parametrize_matrix * point;
        }

        bool is_inside(Eigen::Vector2d point) const { return evaluate_at(point) <= 1; }
};


class EllipseBundle {
    public:
        std::vector<Ellipse> bundle;
        const double x_max, y_max, MW_x, ME_x;
        const double h;

        EllipseBundle(std::map<std::string, double> &geometric_config, const double h_param, const size_t num_ellipses) :
        x_max(geometric_config["x_max"]), y_max(geometric_config["y_max"]), MW_x(geometric_config["MW_x"]),
        ME_x(geometric_config["ME_x"]), h(h_param) {bundle.reserve(num_ellipses);}
            
        void add(const Ellipse &new_ellipse){
            check_intersections(new_ellipse);
            bundle.push_back(new_ellipse);
        }

        const double area() const {
            double total = 0;
            for (auto& ellipse : bundle){
                total += ellipse.area();
            }
            return total;
        }

        Mesh total_mesh() const {
            std::vector<Mesh> mesh_bundle;
            for (const auto &ellipse : bundle) {
                mesh_bundle.push_back(ellipse.get_mesh(h));
            }

            return Mesh::Build(tag::join).meshes(mesh_bundle);
        }

    private:
        // vec1 > vec2 iff at least one of the components is bigger
        bool is_greater(const Eigen::Vector2d &vec1, const Eigen::Vector2d &vec2) const {
            return (vec1.array() > vec2.array()).any();
        }

        void check_intersections(const Ellipse &new_ellipse) const {
            if (!is_inside(new_ellipse)){
                throw std::invalid_argument("Ellipse does not fit in the domain");    
            }
            for (auto & ellipse : bundle){
                if (intersects(new_ellipse, ellipse)){ 
                    throw std::invalid_argument("Ellipses don't have enough gap");
                };
            }
        }
        
        bool is_inside(const Ellipse &e1) const {
            double horizontal_margin = e1.bounds[0] + h;
            double vertical_margin = e1.bounds[1] + h;
            double x = e1.center.x(), y = e1.center.y();

            return ((y > vertical_margin) && (x > horizontal_margin) &&
            (y < y_max - vertical_margin) && (x < x_max - horizontal_margin));
        }

        bool intersects(const Ellipse &e1, const Ellipse &e2) const;
        bool robust_intersect(const Ellipse &e1, const Ellipse &e2) const;
        std::pair<double, double> starting_parameters(const Ellipse &e1, const Ellipse &e2) const;
        std::vector<double> adjacent_angles(const double theta) const;
};

#endif