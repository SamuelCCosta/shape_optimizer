#ifndef ELLIPSE_H
#define ELLIPSE_H

#include "maniFEM.h"
#include "constants.h"
#include "domain_config.h"
#include <numbers>
#include <Eigen/Dense>

class WorkingManifold{
    public:
        WorkingManifold() : current(Manifold::working) {}
        ~WorkingManifold() {current.set_as_working_manifold();}

        WorkingManifold(const WorkingManifold&) = delete;
        WorkingManifold& operator=(const WorkingManifold&) = delete;

        Manifold current;
};


class Ellipse {
    public:
        double height, width, A, B, C, det;
        Eigen::Vector2d center;
        Eigen::Matrix2d M;
        Eigen::Matrix2d transform; //M^-1 = A * A^T <- matriz A

        Ellipse(double x_i, double y_i, double A_i, double B_i, double C_i);

        Mesh get_mesh(const double h, std::list<Manifold>& repository) const;
        Mesh manual_get_mesh(const double h) const;

        double area() const {return pi * 1 / (det * det);}

        Eigen::Vector2d point_at(double theta) const {
            return center + transform * Eigen::Vector2d(std::cos(theta), std::sin(theta));
        }

        Eigen::Vector2d derivative_at(double theta) const {
            return transform * Eigen::Vector2d(-std::sin(theta), std::cos(theta));
        }

        double evaluate_at(const Eigen::Vector2d point) const {
            Eigen::Vector2d d = point - center;
            return d.dot(M * d);
        }
};


class EllipseBundle{
    public:
        std::vector<Ellipse> bundle;
        const DomainConfig cfg;

        EllipseBundle(const DomainConfig& Cfg) : cfg(Cfg) {bundle.reserve(cfg.num_ellipses);}

        void add(const Ellipse &new_ellipse){
            check_intersections(new_ellipse);
            bundle.push_back(new_ellipse);
        }

        void add(Ellipse&& new_ellipse) {
            check_intersections(new_ellipse);
            bundle.push_back(std::move(new_ellipse));
        }

        const double area() const {
            double total = 0;
            for (auto& ellipse : bundle){
                total += ellipse.area();
            }
            return total;
        }

        Mesh total_mesh(std::list<Manifold>& repository) const;
        Mesh manual_total_mesh() const;

    private:
        void check_intersections(const Ellipse &new_ellipse){
            if (!is_inside(new_ellipse)){
                throw std::invalid_argument("Ellipse does not fit in the domain");
            }
            for (auto & ellipse : bundle){
                if (intersects(new_ellipse, ellipse)){ 
                    throw std::invalid_argument("Ellipses don't have enough gap");
                };
            }
        }

        bool is_inside(const Ellipse& e1){ //bounds checking
            double horizontal_margin = cfg.h + e1.width;
            double vertical_margin = cfg.h + e1.height;
            double x = e1.center[0], y = e1.center[1];

            if ((y < vertical_margin) || (x < horizontal_margin) ||
            (y > cfg.y_max - vertical_margin) || (x > cfg.x_max - horizontal_margin) ) {
                return false;
            }
            return true;
        }

        bool intersects(const Ellipse &e1, const Ellipse &e2);
        bool robust_intersect(const Ellipse &e1, const Ellipse &e2) const;
        std::pair<double, double> get_initial_thetas(const Ellipse& e1, const Ellipse& e2) const;
};
#endif