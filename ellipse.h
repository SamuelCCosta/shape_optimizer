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


class Ellipse {
    public:
        double height, width, A, B, C, det;
        Eigen::Vector2d center;
        Eigen::Matrix2d M;
        mutable std::optional<Manifold> representation;
        mutable std::optional<Mesh> mesh;
        

        Ellipse(double x_i, double y_i, double A_i, double B_i, double C_i);

        void meshify() const;

        Mesh get_mesh() const {
            if (! mesh.has_value()){
                meshify();
            }
            return mesh.value();
        }

        double area(){return pi * det;}

        const Eigen::Matrix2d& get_transform_matrix() const {
            if (!transform.has_value()){
                // Cholesky decomposition of M^-1
                Eigen::Matrix2d M_inv = M.inverse();
                Eigen::LLT<Eigen::Matrix2d> llt_inv(M_inv);
                transform = llt_inv.matrixL();
            }
            return transform.value();
        }


        Eigen::Vector2d point_at(double theta) const {
            const auto& L = get_transform_matrix();
            return center + L * Eigen::Vector2d(std::cos(theta), std::sin(theta));
        }

        Eigen::Vector2d derivative_at(double theta) const {
            const auto& L = get_transform_matrix();
            return L * Eigen::Vector2d(-std::sin(theta), std::cos(theta));
        }
        
    private:
        mutable std::optional<Eigen::Matrix2d> transform; //M^-1 = A * A^T <- matriz A
};


class EllipseBundle{
    public:
        std::vector<Ellipse> bundle;
        
        EllipseBundle(){bundle.reserve(num_ellipses);}

        void add(const Ellipse &new_ellipse){
            check_intersections(new_ellipse);
            bundle.push_back(new_ellipse);
        }

        void add_ellipse(Ellipse&& new_ellipse) {
            check_intersections(new_ellipse);
            bundle.push_back(std::move(new_ellipse));
        }

        double area(){
            double total = 0;
            for (auto& ellipse : bundle){
                total += ellipse.area();
            }
            return total;
        }

        Mesh total_mesh();

    private:
        void check_intersections(const Ellipse &new_ellipse){
            for (auto & ellipse : bundle){
                if (intersects(new_ellipse, ellipse)){ 
                    throw std::invalid_argument("Ellipses don't have enough gap");
                };
            }
        };

        bool intersects(const Ellipse &Ellipse1, const Ellipse &Ellipse2);
        //há edge cases em que a verificação simples não dá h distância entre elipses
        bool robust_intersect(const Ellipse &Ellipse1, const Ellipse &Ellipse2) const;
        std::pair<double, double> get_initial_thetas(const Ellipse& e1, const Ellipse& e2) const;
};
#endif