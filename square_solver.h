#ifndef SQUARE_SOLVER_H
#define SQUARE_SOLVER_H

#include "maniFEM.h"
#include "maniUtils.h"
#include "maniSolver.h"
#include "ellipse.h"
#include "constants.h"

#include <cmath>
#include <fstream>
#include <optional>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

using namespace maniFEM;

inline int n_segments(double l){ return std::ceil(l/h); }

double objective(const double heat_sources, const double base_temp, const EllipseBundle& ellipses, bool export_mesh);

//Work in Progress
class SquareSolver {
    public:
        double x_max, y_max, MW_x, ME_x, h;
        double heat_sources, base_temp;
        size_t num_ellipses;
        Mesh sources, south, square_boundary, north;
        const bool export_mesh;
        Manifold ambient;

        SquareSolver(const double x, const double y, const double MW, const double ME, const double heat_sources_,
            const double base_temp_, const double h, const size_t n_ellipses, const bool export_mesh);
        double solve(EllipseBundle& bundle);
};
#endif