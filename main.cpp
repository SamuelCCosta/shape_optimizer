#include "maniFEM.h"
#include "maniUtils.h"
#include "maniSolver.h"
#include "ellipse.h"
#include "constants.h"

#include <cmath>
#include <fstream>
#include <optional>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

using namespace maniFEM;

inline constexpr int n_segments(double l){
    return std::ceil(l/h);
}

int main(){
    bool solve = true;

    Manifold RR2(tag:: Euclid, tag::of_dim, 2);
    Function xy = RR2.build_coordinate_system(tag::Lagrange, tag::of_degree, 1);
    Function x = xy[0], y = xy[1];

    const Function f = 0;
    const Function heat_source = 10; //condição neumann fronteira superior
    const Function base_temp = 0; //condição dirichlet base


    //Quadrado Completo (fronteira)
    const Cell NW(tag::vertex); x(NW) = 0; y(NW) = 1;
    const Cell NE(tag::vertex); x(NE) = 1; y(NE) = 1;
    const Cell SW(tag::vertex); x(SW) = 0; y(SW) = 0;
    const Cell SE(tag::vertex); x(SE) = 1; y(SE) = 0;

    static constexpr float MW_x = 0.3f;
    static constexpr float ME_x = 0.7f;
    const Cell MW(tag::vertex); x(MW) = MW_x; y(MW) = 1;
    const Cell ME(tag::vertex); x(ME) = ME_x; y(ME) = 1;

    constexpr float mid_l = ME_x - MW_x;
    constexpr auto &midW_l = MW_x;
    constexpr float midE_l = 1 - ME_x;
    
    const Mesh MNE = Mesh::Build(tag::grid).shape(tag::segment).start_at(NE).stop_at(ME).divided_in(n_segments(midE_l));
    const Mesh north_middle = Mesh::Build(tag::grid).shape(tag::segment).start_at(ME).stop_at(MW).divided_in(n_segments(mid_l));
    const Mesh NWM = Mesh::Build(tag::grid).shape(tag::segment).start_at(MW).stop_at(NW).divided_in(n_segments(midW_l));
    const Mesh west = Mesh::Build(tag::grid).shape(tag::segment).start_at(NW).stop_at(SW).divided_in(n_segments(1));
    const Mesh south = Mesh::Build(tag::grid).shape(tag::segment).start_at(SW).stop_at(SE).divided_in(n_segments(1));
    const Mesh east = Mesh::Build(tag::grid).shape(tag::segment).start_at(SE).stop_at(NE).divided_in(n_segments(1));

    const Mesh null_neumann = Mesh::Build(tag::join).meshes({east, west, north_middle}); //podemos ignorar
    const Mesh sources = Mesh::Build(tag::join).meshes({NWM,MNE});
    const Mesh square_boundary = Mesh::Build(tag::join).meshes({south, null_neumann, sources});

    //const Mesh square_boundary = Mesh::Build(tag::join).meshes({NWM, north_middle, MNE, east, south, west});

    Ellipse hole(0.4, 0.6, 0.3, 0.1, pi / 3.0);
    
    Mesh hole_mesh = hole.get_mesh();
    
    const Mesh boundary = Mesh::Build(tag::join).mesh(square_boundary).mesh(hole_mesh);

    RR2.set_as_working_manifold(); //não devia ser necessário...

    const Mesh domain = Mesh::Build(tag::frontal).boundary(boundary).desired_length(h);

    if (!solve) {domain.export_to_file(tag::gmsh, "domain.msh");}

    if (solve){
    std::map<Cell, size_t> numbering = create_node_numbering(domain, degree);
    
    Eigen::VectorXd solution = (hand_coded)? build_poisson_solution(f, base_temp, heat_source, domain, south, sources, numbering, tag::hand_coded):
                                             build_poisson_solution(f, base_temp, heat_source, domain, south, sources, square_boundary, numbering, degree);
    
    if (degree == 1 && export_file) {
        domain.export_to_file ( tag::gmsh, "solution.msh", numbering );

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
        solution_file << domain.number_of ( tag::vertices ) << std::endl;  // number of values listed below
        Mesh::Iterator it = domain.iterator ( tag::over_vertices );
        for ( it .reset(); it .in_range(); it++ )
        {	Cell P = *it;
            const size_t i = numbering [P];
            solution_file << i+1 << " " << solution [i] << std::endl;   }
        }
    }
    
    } //if solve

	return 0;
}

// scp samuel@192.168.1.145:/home/samuel/shape_optimizer/domain.msh .