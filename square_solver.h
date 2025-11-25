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

inline constexpr int n_segments(double l){ return std::ceil(l/h); }

double objective(const double heat_sources, const double base_temp, const EllipseBundle& ellipses, bool export_mesh);

//Work in Progress
/*
class SquareSolver {
    public:
        double x_max, y_max, MW_x, ME_x, h;
        size_t num_ellipses;

        SquareSolver(double x, double y, double MW, double ME, double h, size_t n_ellipses);
        void solve();
}; */
#endif