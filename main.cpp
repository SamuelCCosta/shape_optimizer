#include "square_solver.h"
#include <chrono>

using namespace maniFEM;

void objective_ellipses(){
    //Manifold RR2(tag:: Euclid, tag::of_dim, 2);
    //Function xy = RR2.build_coordinate_system(tag::Lagrange, tag::of_degree, 1);
    //Function x = xy[0], y = xy[1];

    const double heat_source = 10; //condição neumann fronteira superior
    const double base_temp = 0; //condição dirichlet na base
    constexpr bool export_mesh = false;

    EllipseBundle ellipses;

    ellipses.add(Ellipse(0.5, 0.5, 78.1888, -37.7726, 34.5663)); // a = 0.28, b = 0.1, theta = pi/3
    ellipses.add(Ellipse(0.7, 0.31, 204.082, 0.0, 12.755)); // a = 0.28, b = 0.07, theta = pi/2
    ellipses.add(Ellipse(0.19, 0.59, 45.0817, -25.1834, 200.0929)); //a = 0.156, b = 0.07, theta = pi/20
    ellipses.add(Ellipse(0.19, 0.17, 45.0817, -25.1834, 200.0929)); //a = 0.156, b = 0.07, theta = pi/20

    double value = objective(heat_source, base_temp, ellipses, export_mesh);
    std::cout << "Objective: " << value << std::endl;
}

int main(){
    auto start = std::chrono::high_resolution_clock::now();

    objective_ellipses();

    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;
    return 0;
}

// para h = 0.02, o código é executado em ~170ms (Ryzen 5 7600)
// sem escrever os ficheiros .msh, demora ~100ms