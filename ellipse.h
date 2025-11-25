#ifndef ELLIPSE_H
#define ELLIPSE_H

#include "maniFEM.h"
#include "constants.h"
#include <optional>
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


class alignas(256) Ellipse {
    public:
        double height, width, A, B, C, det;
        Eigen::Vector2d center;
        Eigen::Matrix2d M;
        Eigen::Matrix2d transform; //M^-1 = A * A^T <- matriz A
        mutable Manifold representation;
        mutable Mesh mesh;

        Ellipse(double x_i, double y_i, double A_i, double B_i, double C_i);

        void meshify() const;

        Mesh get_mesh() const {
            if (!mesh.exists()){
                meshify();
            }
            return mesh;
        }

        const double area() const {return pi * 1/ (det * det);}

        Eigen::Vector2d point_at(double theta) const {
            return center + transform * Eigen::Vector2d(std::cos(theta), std::sin(theta));
        }

        Eigen::Vector2d derivative_at(double theta) const {
            return transform * Eigen::Vector2d(-std::sin(theta), std::cos(theta));
        }
};


class EllipseBundle{
    public:
        std::vector<Ellipse> bundle;
        
        EllipseBundle(){bundle.reserve(num_ellipses);}

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

        const Mesh total_mesh() const;

    private:
        void check_intersections(const Ellipse &new_ellipse){
            for (auto & ellipse : bundle){
                if (intersects(new_ellipse, ellipse)){ 
                    throw std::invalid_argument("Ellipses don't have enough gap");
                };
            }
        };

        bool intersects(const Ellipse &e1, const Ellipse &e2);
        bool robust_intersect(const Ellipse &e1, const Ellipse &e2) const;
        std::pair<double, double> get_initial_thetas(const Ellipse& e1, const Ellipse& e2) const;
};
#endif