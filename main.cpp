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

    SquareSolver sqs = SquareSolver(geometric_info, h, heat_source, base_temp, true, false);

    EllipseBundle bundle(geometric_info, h, num_ellipses);
    /*
    bundle.add(Ellipse(0.8618206759053694, 0.638244200440503, 225.67493248788125, -160.57422243346917, 280.5987404955669));
    bundle.add(Ellipse(0.8379835207017033, 0.3939957829236339, 192.45436487876822, -21.333606792881376, 78.24251380072735));
    bundle.add(Ellipse(0.15718368156532972, 0.11260201636280878, 398.10059552791836, -241.77576084751604, 451.12919820660574));
    bundle.add(Ellipse(0.3854623436369899, 0.4835325070889105, 14.240792213872059, 7.025599246051849, 8.336150579822723));
    */
    //bundle.add(Ellipse(0.3, 0.3, 81.0, 15.0, 81.0));

    //bundle.add(Ellipse(0.33, 0.33, 11.111111, 0.0, 69.444444));
    //bundle.add(Ellipse(0.675, 0.675, 69.444444, 0.0, 11.111111));

    bundle.add(Ellipse(0.6761914390629538, 0.8331474909659954, 234.62634533093578, 55.20225964246275, 294.41739831068963));
    bundle.add(Ellipse(0.14177601073689275, 0.38267316545502855, 166.83879221822332, -70.30573735061755, 214.1889264220846));
    bundle.add(Ellipse(0.38952921827986187, 0.4131190477291755, 178.83633814702222, -3.5752206604725445, 70.3255660637861));
    bundle.add(Ellipse(0.8981855297026812, 0.696113885664253, 241.26170879026645, -26.515173296621807, 269.63343984075857));
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
