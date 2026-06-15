#include "square_solver.h"
#include <chrono>
#include <iostream>

using namespace maniFEM;

void objective_ellipses(){
    const double h = 0.001;
    std::cout << "h: " << h << std::endl;
    const double heat_source = 10.0; //condição neumann fronteira superior
    const double base_temp = 0.0; //condição dirichlet na base
    const size_t num_ellipses = 4; //max number of ellipses

    std::map<std::string, double> geometric_info = {{"x_max", 1.0}, {"y_max", 1.0}, {"MW_x", 0.3}, {"ME_x", 0.7}};

    //bool export_domain = false, export_solution = false;
    //SquareSolver sqs = SquareSolver(geometric_info, h, heat_source, base_temp, export_domain, export_solution);

    SquareSolver sqs = SquareSolver(geometric_info, h, heat_source, base_temp, false, false);

    EllipseBundle bundle(geometric_info, h, num_ellipses);
    
    //bundle.generate_random(1234);

    const double pen = 0.0;
    double area_percent = 1-bundle.area();
    bool print_area = true;
    if (print_area) { std::cout << "Total area: " << area_percent << std::endl; }
    double final_pen = pen * area_percent;

    double final_result = sqs.solve(bundle);

    std::cout << "No pen: " << final_result << std::endl;
    if (pen != 0.0) {
        std::cout << "Penalization: " << final_pen << std::endl;
        std::cout << "Penalized: " << final_result + final_pen << std::endl;
    }
    if (area_percent == 1.0) { 
        std::cout << "Diff between real and calculated: " << 3.9174195203606 - final_result << std::endl;
    }
}

int main(){
    auto start = std::chrono::high_resolution_clock::now();

    objective_ellipses();    

    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;
    return 0;
}
