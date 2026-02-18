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
    
    for(int i = 0; i < 1; i++){
        EllipseBundle ellipses(cfg);


        double value = sqs.solve(ellipses);
        std::cout << "Objective: " << value << std::endl;         
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

// para h = 0.02, o código é executado em ~170ms (Ryzen 5 7600)
// sem escrever os ficheiros .msh, demora ~100ms
