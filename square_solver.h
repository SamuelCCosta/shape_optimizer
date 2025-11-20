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

double objective(const Function heat_sources, const Function base_temp, const EllipseBundle& ellipses, bool export_mesh);
#endif