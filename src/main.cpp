#include "square_solver.h"
#include <chrono>
#include <iostream>
#include <limits>
#include <iomanip>

using namespace maniFEM;

double objective_ellipses(const double h){
    constexpr bool COUT_OBJECTIVE_ELLIPSES = false;
    if (COUT_OBJECTIVE_ELLIPSES) { std::cout << "h: " << h << std::endl; }
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
    if (print_area && COUT_OBJECTIVE_ELLIPSES) { std::cout << "Total area: " << area_percent << std::endl; }
    double final_pen = pen * area_percent;

    double final_result = sqs.solve(bundle);

    if(COUT_OBJECTIVE_ELLIPSES){ std::cout << "No pen: " << final_result << std::endl; }
    
    if (pen != 0.0 && COUT_OBJECTIVE_ELLIPSES) {
        std::cout << "Penalization: " << final_pen << std::endl;
        std::cout << "Penalized: " << final_result + final_pen << std::endl;
    }
    if (area_percent == 1.0 && COUT_OBJECTIVE_ELLIPSES) { 
        std::cout << "Diff between real and calculated: " << 3.917419523584066 - final_result << std::endl;
    }
    return final_result;
}

int main(){
    const double h_max = 0.1, h_min = 0.001;
    const int h_n = 15;
    std::vector<double> h_tests;
    for (int k = 0; k < h_n; k++) {
        double h_k = h_max * std::pow(h_min/h_max, static_cast<double>(k) / (h_n-1));
        h_tests.push_back(h_k); 
    }

    for (const auto &h_k : h_tests) {
        std::cout << "h: " << h_k;
        auto start = std::chrono::high_resolution_clock::now();

        double final_result = objective_ellipses(h_k);

        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        //std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;
        std::cout << std::setprecision(std::numeric_limits<double>::max_digits10)
          << " val: " << final_result
          << " exec time: " << duration
          << std::endl;
    }
    return 0;
}
