#include "square_solver.h"

/*
Manifold fallback(tag::non_existent);
    if (!Manifold::working.exists()) {
        fallback = Manifold(tag::Euclid, tag::of_dim, 2);
        Manifold::working.build_coordinate_system(tag::Lagrange, tag::of_degree, 1);
    }
    Manifold& active_manifold = Manifold::working;
    Function xy = active_manifold.coordinates();
    Function x = xy[0], y = xy[1];
*/

SquareSolver::SquareSolver(DomainConfig Cfg, const double heat_sources_, const double base_temp_, const double penalization_, const bool export_mesh_) :
    cfg(Cfg), heat_sources(heat_sources_), base_temp(base_temp_), penalization(penalization_), sources(tag::non_existent), south(tag::non_existent),
    square_boundary(tag::non_existent), north(tag::non_existent), export_mesh(export_mesh_), ambient(tag::non_existent){
    ambient = Manifold(tag::Euclid, tag::of_dim, 2);
    Function xy = ambient.build_coordinate_system(tag::Lagrange, tag::of_degree, 1);
    Function x = xy[0], y = xy[1];

    //Quadrado Completo (fronteira)
    const Cell NW(tag::vertex); x(NW) = 0; y(NW) = cfg.y_max;
    const Cell NE(tag::vertex); x(NE) = cfg.x_max; y(NE) = cfg.y_max;
    const Cell SW(tag::vertex); x(SW) = 0; y(SW) = 0;
    const Cell SE(tag::vertex); x(SE) = cfg.x_max; y(SE) = 0;

    const Cell MW(tag::vertex); x(MW) = cfg.MW_x; y(MW) = cfg.y_max;
    const Cell ME(tag::vertex); x(ME) = cfg.ME_x; y(ME) = cfg.y_max;

    double mid_l = cfg.ME_x - cfg.MW_x;
    auto &midW_l = cfg.MW_x;
    double midE_l = cfg.x_max - cfg.ME_x;
    
    const Mesh MNE = Mesh::Build(tag::grid).shape(tag::segment).start_at(NE).stop_at(ME).divided_in(cfg.n_segments(midE_l));
    const Mesh north_middle = Mesh::Build(tag::grid).shape(tag::segment).start_at(ME).stop_at(MW).divided_in(cfg.n_segments(mid_l));
    const Mesh NWM = Mesh::Build(tag::grid).shape(tag::segment).start_at(MW).stop_at(NW).divided_in(cfg.n_segments(midW_l));
    const Mesh west = Mesh::Build(tag::grid).shape(tag::segment).start_at(NW).stop_at(SW).divided_in(cfg.n_segments(cfg.y_max));
    south = Mesh::Build(tag::grid).shape(tag::segment).start_at(SW).stop_at(SE).divided_in(cfg.n_segments(cfg.x_max));
    const Mesh east = Mesh::Build(tag::grid).shape(tag::segment).start_at(SE).stop_at(NE).divided_in(cfg.n_segments(cfg.y_max));

    const Mesh null_neumann = Mesh::Build(tag::join).meshes({east, west, north_middle}); //podemos ignorar na construção das condições
    sources = Mesh::Build(tag::join).meshes({NWM,MNE});
    square_boundary = Mesh::Build(tag::join).meshes({south, null_neumann, sources});
    north = Mesh::Build(tag::join).meshes({sources, north_middle});
}

double SquareSolver::solve(EllipseBundle& bundle){
    Mesh inner_boundary = bundle.total_mesh();
    
    Mesh boundary = Mesh::Build(tag::join).mesh(square_boundary).mesh(inner_boundary);

    const Mesh domain = Mesh::Build(tag::frontal).boundary(boundary).desired_length(cfg.h);

    std::map<Cell, size_t> numbering = create_node_numbering(domain, degree);
    
    Eigen::VectorXd solution = build_laplace_solution(base_temp, heat_sources, domain, south, sources, numbering);

    //função objetivo
    double objective = 0;
    double volume = cfg.x_max * cfg.y_max - bundle.area();

    //Podemos implementar um loop simples mais tarde, integrar trapézios é fácil
    // J(u) = integrate u in north + penalization *volume
    {
    FiniteElement fe_bdry(tag::with_master, tag::segment, tag::Lagrange, tag::of_degree, 1);
    Integrator integr_bdry = fe_bdry.set_integrator(tag::Gauss, tag::seg_3);
    
    Mesh::Iterator it = north.iterator(tag::over_cells_of_max_dim);
    for(it.reset(); it.in_range(); it++){
        Cell seg = *it;
        fe_bdry.dock_on(seg);
        Mesh::Iterator it_vert = seg.boundary().iterator(tag::over_vertices, tag::force_positive);
        for(it_vert.reset(); it_vert.in_range(); it_vert++){
            Cell V = *it_vert;
            Function phi_V = fe_bdry.basis_function(V);
            objective +=
                fe_bdry.integrate(solution[numbering.at(V)] * phi_V);
        }
    }
    }

    double final_objective = objective + penalization * volume;

    if (export_mesh){
        domain.export_to_file(tag::gmsh, "domain.msh");
        
        domain.export_to_file (tag::gmsh, "solution.msh", numbering);

        //Incluir dados da solução no .msh
        {
        std::ofstream solution_file ("solution.msh", std::fstream::app);
        solution_file << "$NodeData" << std::endl;
        solution_file << "1" << std::endl;   // one string follows
        solution_file << "\"Solution\"" << std::endl;
        solution_file << "1" << std::endl;   //  one real follows
        solution_file << "0.0" << std::endl;  // time [??]
        solution_file << "3" << std::endl;   // three integers follow
        solution_file << "0" << std::endl;   // time step [??]
        solution_file << "1" << std::endl;  // scalar values of u
        solution_file << domain.number_of (tag::vertices) << std::endl;  // number of values listed below
        Mesh::Iterator it = domain.iterator (tag::over_vertices);
        for (it .reset(); it .in_range(); it++)
        {	Cell P = *it;
            const size_t i = numbering [P];
            solution_file << i+1 << " " << solution [i] << std::endl;   }
        }
    }

    return final_objective;
}