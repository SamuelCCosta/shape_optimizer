#include "square_solver.h"
#include <chrono>
#include <iostream>

using namespace maniFEM;

void objective_ellipses(){
    const double h = 0.02;

    const double heat_source = 10.0; //condição neumann fronteira superior
    const double base_temp = 0.0; //condição dirichlet na base
    const size_t num_ellipses = 4; //max number of ellipses

    std::map<std::string, double> geometric_info = {{"x_max", 1.0}, {"y_max", 1.0}, {"MW_x", 0.3}, {"ME_x", 0.7}};

    //bool export_domain = false, export_solution = false;
    //SquareSolver sqs = SquareSolver(geometric_info, h, heat_source, base_temp, export_domain, export_solution);

    SquareSolver sqs = SquareSolver(geometric_info, h, heat_source, base_temp, true, true);

    EllipseBundle bundle(geometric_info, h, num_ellipses);
    
    bundle.add(Ellipse(0.3, 0.5, 27.777, -16.666, 27.777));
    bundle.add(Ellipse(0.7, 0.5, 27.777, -16.666, 27.777));

    /*bundle.add(Ellipse(0.53446306, 0.50259374, 243.77396884, -33.96889496, 261.82802112));
    bundle.add(Ellipse(0.42235164, 0.70059529, 233.99220896, -43.95489961, 288.09370986));
    bundle.add(Ellipse(0.56719944, 0.34218724, 289.65241533, -106.94255163, 203.10652993));
    bundle.add(Ellipse(0.46513083, 0.8759024, 208.72228596, 71.36566829, 130.41681181));*/

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
