#ifndef SQUARE_SOLVER_H
#define SQUARE_SOLVER_H

#include "maniFEM.h"
#include "ellipse.h"

#include <cmath>
#include <fstream>
#include <optional>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

//temp
#include <chrono>

using namespace maniFEM;


class SquareSolver {
    public:
        double x_max, y_max, MW_x, ME_x;
        double h;
        double heat_sources, base_temp; //boundary conditions
        const bool export_domain;
        const bool export_result;
        Manifold ambient{tag::Euclid, tag::of_dim, 2}; //R^2
        Mesh north{tag::non_existent};
        Mesh south{tag::non_existent}; //Dirichlet boundary condition
        Mesh sources{tag::non_existent}; //Neumann boundary condition and objective funcional calculation
        Mesh square_boundary{tag::non_existent}; //Domain construction

        SquareSolver(std::map<std::string, double> &geometric_config, const double h_param, const double heat_srcs,
            const double base_tmp, const bool export_dom, const bool export_res);
        
        // Both exports are set to false
        SquareSolver(std::map<std::string, double> &geometric_config, const double h_param, const double heat_srcs,
            const double base_tmp);

        // returns the objective functional value
        double solve(EllipseBundle &bundle);

        double solve_frontal(EllipseBundle &bundle);

    private:
        std::map<Cell, size_t> create_numbering(const Mesh &domain);

        //specialization for constant functions
        Eigen::VectorXd build_laplace_solution(const Mesh &domain, const std::map<Cell,size_t> &numbering);

        void impose_value_of_unknown(Eigen::SparseMatrix<double> &matrix_A, 
            Eigen::VectorXd &vector_b, const size_t cell_id, const double val);

        double objective_no_penalty(const Eigen::VectorXd &solution, const std::map<Cell,size_t> &numbering);

        int n_segments(double l){ return std::ceil(l/h); }

        int n_segments(double l, double size){ return std::ceil(l/size); }

        std::vector<Eigen::Vector2d> get_delaunay_grid(EllipseBundle &bundle);
};

#endif