#include "maniFEM.h"
#include "maniUtils.h"
#include "maniSolver.h"

#include <fstream>
#include <optional>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

using namespace maniFEM;

int main(){
    // Constantes
    const double h = 0.1; // comprimento médio da malha
    const double radius = 1.; // raio do disco
    constexpr size_t degree = 1;
    constexpr bool hand_coded = 1;
    const bool export_file = 1;
    assert(degree == 1 || degree == 2);
    assert(!(degree == 2 && hand_coded));


    Manifold RR2(tag:: Euclid, tag::of_dim, 2);
    Function xy = RR2.build_coordinate_system(tag::Lagrange, tag::of_degree, 1);
    Function x = xy[0], y = xy[1];


    /* u(x,y) = x^2 + y^2  
    -lap u = -4 u(S^1)= 1 du/dn(S1) = 2
    f = -4 , g = 1, normal_deriv = 2
    */

    //const Function f = -1 * exp(-x*x-y*y) * sin(x*y) * cos(x*y); // -lap u = f
    const Function f = -4;

    const Function g = 1; //Condição de fronteira de Dirichlet

    const Function normal_deriv = 2; //condição de fronteira de Neumann

    const Manifold circle = RR2.implicit(x*x + y*y == radius * radius);

    //Fronteira Dirichlet/Neumann
    const Cell start_dirichlet(tag::vertex); x(start_dirichlet) = radius; y(start_dirichlet) = 0;
    const Cell end_dirichlet(tag::vertex); x(end_dirichlet) = -radius; y(end_dirichlet) = 0;
    const Cell& start_neumann = end_dirichlet; 
    const Cell& end_neumann = start_dirichlet;
    const std::vector<double> dirichlet_direction = {0,1};
    const std::vector<double> neumann_direction = {0,-1};

    Mesh dirichlet_bdry = Mesh::Build(tag::frontal).start_at(start_dirichlet)
                        .towards(dirichlet_direction).stop_at(end_dirichlet).desired_length(h);
    Mesh neumann_bdry = Mesh::Build(tag::frontal).start_at(start_neumann)
                            .towards(neumann_direction).stop_at(end_neumann).desired_length(h);

    Mesh bdry = Mesh::Build(tag::join).meshes({dirichlet_bdry,neumann_bdry});

    RR2.set_as_working_manifold();
    Mesh disk = Mesh::Build (tag::frontal). boundary(bdry). desired_length(h);
    
    std::map<Cell, size_t> numbering = create_node_numbering(disk, degree);
    
    
    Eigen::VectorXd solution = (hand_coded)? build_poisson_solution(f, g, normal_deriv, disk, dirichlet_bdry, neumann_bdry, numbering, tag::hand_coded):
                                             build_poisson_solution(f, g, normal_deriv, disk, dirichlet_bdry, neumann_bdry, bdry, numbering, degree);
    

    //Exportar ficheiro
    if (degree == 1 && export_file) {
        disk.export_to_file ( tag::gmsh, "solution.msh", numbering );

        //Mexer no .msh
        {
        std::ofstream solution_file ( "solution.msh", std::fstream::app );
        solution_file << "$NodeData" << std::endl;
        solution_file << "1" << std::endl;   // one string follows
        solution_file << "\"Solution\"" << std::endl;
        solution_file << "1" << std::endl;   //  one real follows
        solution_file << "0.0" << std::endl;  // time [??]
        solution_file << "3" << std::endl;   // three integers follow
        solution_file << "0" << std::endl;   // time step [??]
        solution_file << "1" << std::endl;  // scalar values of u
        solution_file << disk.number_of ( tag::vertices ) << std::endl;  // number of values listed below
        Mesh::Iterator it = disk.iterator ( tag::over_vertices );
        for ( it .reset(); it .in_range(); it++ )
        {	Cell P = *it;
            const size_t i = numbering [P];
            solution_file << i+1 << " " << solution [i] << std::endl;   }
        }
    }
    
	return 0;
}