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
        
        ellipses.add(Ellipse(0.5, 0.5, 78.1888, -37.7726, 34.5663)); // a = 0.28, b = 0.1, theta = pi/3
        ellipses.add(Ellipse(0.7, 0.31, 204.082, 0.0, 12.755)); // a = 0.28, b = 0.07, theta = pi/2
        ellipses.add(Ellipse(0.19, 0.59, 45.0817, -25.1834, 200.0929)); //a = 0.156, b = 0.07, theta = pi/20
        ellipses.add(Ellipse(0.19, 0.17, 45.0817, -25.1834, 200.0929)); //a = 0.156, b = 0.07, theta = pi/20
        
        double value = sqs.solve(ellipses);
        std::cout << "Objective: " << value << std::endl;
        
        /*
        EllipseBundle ellipses2(cfg);
        ellipses2.add(Ellipse(0.09846759459531274, 0.8779145670536423, 360.2709520038292, 112.38584975360911, 393.63233210725235));
        ellipses2.add(Ellipse(0.5056591695680825, 0.3407869436316506, 162.84307818263258, -62.60305881718613, 274.44694581695825));
        ellipses2.add(Ellipse(0.766743884499224, 0.46129767698468077, 186.3704649771612, -144.589190636172, 227.2284762469888));
        ellipses2.add(Ellipse(0.7819328169643374, 0.6512709822028184, 309.4108476099627, 195.55068464492143, 338.3886609589487));

        double value2 = sqs.solve(ellipses2);
        std::cout << "Objective: " << value2 << std::endl;
        */
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


/* O FiniteElement não quer funcionar (ver malha em problematic_domain.msh)
EllipseBundle ellipses2(cfg);
ellipses2.add(Ellipse(0.6233130295388306, 0.697342662933589, 174.0950667024922, 159.36080370536146, 180.09440423275473));
ellipses2.add(Ellipse(0.15383334128618065, 0.3584060814932851, 302.79239251802045, -65.13477823510976, 277.9225236672475));
ellipses2.add(Ellipse(0.45328057400325505, 0.4218320070826215, 149.41014091251955, -80.78375075477966, 325.7908186272372));
ellipses2.add(Ellipse(0.6902436338368403, 0.41948604161280967, 337.44222727593245, -145.11065908780162, 209.55431092846715));
*/