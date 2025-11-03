#ifndef FEM_UTILS
#define FEM_UTILS

#include "maniFEM.h"

#include <utility> // For std::pair
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

using namespace maniFEM;


void impose_value_of_unknown(Eigen::SparseMatrix<double>& matrix_A, 
    Eigen::VectorXd& vector_b, const size_t i, const double val);

std::pair<Mesh, Mesh> create_mesh_from_implicit(Manifold& ambient_manifold, 
    const Function& implicit_eq, const double h, const Cell& start_point);

std::map<Cell, size_t> create_node_numbering(const Mesh& mesh, const size_t degree);

#endif