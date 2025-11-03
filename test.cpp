#include "maniFEM.h"
#include "maniUtils.h"

#include <fstream>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

using namespace maniFEM;


int main(){
    // Constantes
    const double h = 0.1; // comprimento médio da malha
    const double radius = 1.; // raio do disco

    // criar R2 e definir coordenadas
    Manifold RR2(tag:: Euclid, tag::of_dim, 2);
    Function xy = RR2.build_coordinate_system(tag::Lagrange, tag::of_degree, 1);
    Function x = xy[0], y = xy[1];

    //Implicit e ponto de partida da malha
    const Manifold circle = RR2.implicit(x*x + y*y == radius * radius); //conjunto de nível 0

    /*Neste exemplo, a condição de dirichlet é na parte de cima da circunferência e a de neumann em baixo*/
    const Cell start_dirichlet(tag::vertex); x(start_dirichlet) = radius; y(start_dirichlet) = 0; //o fim do neumann
    const Cell end_dirichlet(tag::vertex); x(end_dirichlet) = -radius; y(end_dirichlet) = 0; //também o início do neumann
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

    Mesh::Composite domain = Mesh::Build(tag::gather).cell("start dirichlet", start_dirichlet)
                .cell("end dirichlet", end_dirichlet).mesh("dirichlet", dirichlet_bdry)
                .mesh("neumann", neumann_bdry).mesh("boundary", bdry).mesh("disk", disk);

    domain.export_to_file(tag::gmsh, "domain.msh");
    return 0;
}