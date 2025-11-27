#include "square_solver.h"
#include <chrono>

using namespace maniFEM;

void objective_ellipses(){
    const double h = 0.02;
    const double heat_source = 10; //condição neumann fronteira superior
    const double base_temp = 0; //condição dirichlet na base
    const double penalization = 10; //penalização no volume

    bool export_mesh = false;
    DomainConfig cfg = DomainConfig(1, 1, 0.3, 0.7, h, 10);
    SquareSolver sqs = SquareSolver(cfg, heat_source, base_temp, penalization, export_mesh);
    
    EllipseBundle ellipses(cfg);
    /*
    ellipses.add(Ellipse(0.5, 0.5, 78.1888, -37.7726, 34.5663)); // a = 0.28, b = 0.1, theta = pi/3
    ellipses.add(Ellipse(0.7, 0.31, 204.082, 0.0, 12.755)); // a = 0.28, b = 0.07, theta = pi/2
    ellipses.add(Ellipse(0.19, 0.59, 45.0817, -25.1834, 200.0929)); //a = 0.156, b = 0.07, theta = pi/20
    ellipses.add(Ellipse(0.19, 0.17, 45.0817, -25.1834, 200.0929)); //a = 0.156, b = 0.07, theta = pi/20
    */
    ellipses.add(Ellipse(5.91629030e-01,  6.40340026e-01,  2.04159081e+02, -4.35976084e+01, 1.12926963e+02));
    ellipses.add(Ellipse(5.42716448e-01,  7.59968556e-01,  1.00550021e+02, -1.30988244e+01, 1.17660615e+02));

    //std::cout << ellipses.area() << std::endl;
    //double value = sqs.solve(ellipses);
    //std::cout << "Objective: " << value << std::endl;
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