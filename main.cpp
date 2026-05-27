#include "square_solver.h"
#include <chrono>
#include <iostream>

using namespace maniFEM;

void objective_ellipses(){
    const double h = 0.015;

    const double heat_source = 10.0; //condição neumann fronteira superior
    const double base_temp = 0.0; //condição dirichlet na base
    const size_t num_ellipses = 4; //max number of ellipses

    std::map<std::string, double> geometric_info = {{"x_max", 1.0}, {"y_max", 1.0}, {"MW_x", 0.3}, {"ME_x", 0.7}};

    //bool export_domain = false, export_solution = false;
    //SquareSolver sqs = SquareSolver(geometric_info, h, heat_source, base_temp, export_domain, export_solution);

    SquareSolver sqs = SquareSolver(geometric_info, h, heat_source, base_temp, false, false);

    EllipseBundle bundle(geometric_info, h, num_ellipses);
    /*
    bundle.add(Ellipse(0.8618206759053694, 0.638244200440503, 225.67493248788125, -160.57422243346917, 280.5987404955669));
    bundle.add(Ellipse(0.8379835207017033, 0.3939957829236339, 192.45436487876822, -21.333606792881376, 78.24251380072735));
    bundle.add(Ellipse(0.15718368156532972, 0.11260201636280878, 398.10059552791836, -241.77576084751604, 451.12919820660574));
    bundle.add(Ellipse(0.3854623436369899, 0.4835325070889105, 14.240792213872059, 7.025599246051849, 8.336150579822723));
    */
    //bundle.add(Ellipse(0.3, 0.3, 81.0, 15.0, 81.0));

    const double pen = 0.0;
    double area_percent = 1-bundle.area();
    double final_pen = pen * area_percent;

    double final_result = sqs.solve(bundle);

    std::cout << "No pen: " << final_result << std::endl;
    std::cout << "Penalization: " << final_pen << std::endl;
    std::cout << "Penalized: " << final_result + final_pen << std::endl;
}

int main(){
    auto start = std::chrono::high_resolution_clock::now();

    objective_ellipses();    

    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;
    return 0;
}
