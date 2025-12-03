#ifndef SQUARE_SOLVER_H
#define SQUARE_SOLVER_H

#include "maniFEM.h"
#include "maniUtils.h"
#include "maniSolver.h"
#include "ellipse.h"
#include "constants.h"
#include "domain_config.h"

#include <cmath>
#include <fstream>
#include <optional>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

using namespace maniFEM;

class SquareSolver {
    public:
        DomainConfig cfg;
        double heat_sources, base_temp;
        double penalization;
        Mesh sources, south, square_boundary, north;
        const bool export_mesh;
        Manifold ambient;

        SquareSolver(DomainConfig cfg, const double heat_sources, const double base_temp, const double penalization, const bool export_mesh);
        double solve(EllipseBundle& bundle);
};
#endif