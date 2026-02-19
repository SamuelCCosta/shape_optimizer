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
    bundle.add(Ellipse(0.4528497255555963, 0.5195018072130883, 255.24891698211133, -38.35196654400833, 262.5530482957416));
    //std::cout << "ellipse 1 added" << std::endl;
    bundle.add(Ellipse(0.45017099450198184, 0.6923857353204883, 231.7094153654055, -47.15190162486053, 287.10515606531953));
    //std::cout << "ellipse 2 added" << std::endl;
    bundle.add(Ellipse(0.5265268828586713, 0.33421315950352537, 293.72495839816185, -112.08640183253189, 204.51894409046088));
    //std::cout << "ellipse 3 added" << std::endl;
    bundle.add(Ellipse(0.4319381994260433, 0.8739065566217921, 209.82444979164484, 67.87935546812955, 128.8843669760333));
    //std::cout << "ellipse 4 added" << std::endl;

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
