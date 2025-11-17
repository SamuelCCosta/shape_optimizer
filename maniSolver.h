#ifndef FEM_SOLVER
#define FEM_SOLVER

#include "maniFEM.h"
#include <optional>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

using namespace maniFEM;

Eigen::VectorXd build_poisson_solution(const Function& f, const Function& boundary_value, const Function& normal_deriv,
    const Mesh& domain, const Mesh& dirichlet_boundary, const Mesh& neumann_boundary, const Mesh& boundary, const std::map<Cell, size_t>& numbering, const size_t& degree);

Eigen::VectorXd build_poisson_solution(const Function& f, const Function& boundary_value, const Function& normal_deriv,
    const Mesh& domain, const Mesh& dirichlet_boundary, const Mesh& neumann_boundary, const std::map<Cell, size_t>& numbering, const tag::HandCoded&);

//builds laplace equation solution using hand coded Lagrange P1 finite elements
Eigen::VectorXd build_laplace_solution(const Function& boundary_value, const Function& normal_deriv, const Mesh& domain,
    const Mesh& dirichlet_boundary, const Mesh& neumann_boundary, const std::map<Cell, size_t>& numbering);

Eigen::VectorXd build_laplace_solution(const double boundary_value, const double normal_deriv, const Mesh& domain,
    const Mesh& dirichlet_boundary, const Mesh& neumann_boundary, const std::map<Cell, size_t>& numbering);
#endif