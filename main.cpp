#include "square_solver.h"
#include <chrono>
#include <iostream>

using namespace maniFEM;

void objective_ellipses(){
    const double h = 0.02;
    const double heat_source = 10; //condição neumann fronteira superior
    const double base_temp = 0; //condição dirichlet na base
    const double penalization = 100; //penalização no volume
    const size_t num_ellipses = 5;

    std::map<std::string, double> geometric_info = {{"x_max", 1.0}, {"y_max", 1.0}, {"MW_x", 0.3}, {"ME_x", 0.7}};

    bool export_domain = false, export_solution = false;
    
    SquareSolver sqs = SquareSolver(geometric_info, h, heat_source, base_temp, penalization, export_domain, export_solution);

    EllipseBundle bundle(geometric_info, h, num_ellipses);
    bundle.add(Ellipse(0.5,0.5,13,0,13));

    double final_result = sqs.solve(bundle);

    std::cout << final_result << std::endl;
}

int main(){
    auto start = std::chrono::high_resolution_clock::now();

    objective_ellipses();    

    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;
    return 0;
}
