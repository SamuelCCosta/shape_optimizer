#include "square_solver.h"

//apenas para o quadrado [0,1]
double objective(const double heat_sources, const double base_temp, const EllipseBundle& ellipses, const bool export_mesh){
    Manifold fallback(tag::non_existent);
    if (!Manifold::working.exists()) {
        fallback = Manifold(tag::Euclid, tag::of_dim, 2);
        Manifold::working.build_coordinate_system(tag::Lagrange, tag::of_degree, 1);
    }
    Manifold& active_manifold = Manifold::working;
    Function xy = active_manifold.coordinates();
    Function x = xy[0], y = xy[1];

    //Quadrado Completo (fronteira)
    const Cell NW(tag::vertex); x(NW) = 0; y(NW) = 1;
    const Cell NE(tag::vertex); x(NE) = 1; y(NE) = 1;
    const Cell SW(tag::vertex); x(SW) = 0; y(SW) = 0;
    const Cell SE(tag::vertex); x(SE) = 1; y(SE) = 0;

    static constexpr double MW_x = 0.3f;
    static constexpr double ME_x = 0.7f;
    const Cell MW(tag::vertex); x(MW) = MW_x; y(MW) = 1;
    const Cell ME(tag::vertex); x(ME) = ME_x; y(ME) = 1;

    constexpr double mid_l = ME_x - MW_x;
    constexpr auto &midW_l = MW_x;
    constexpr double midE_l = 1 - ME_x;
    
    const Mesh MNE = Mesh::Build(tag::grid).shape(tag::segment).start_at(NE).stop_at(ME).divided_in(n_segments(midE_l));
    const Mesh north_middle = Mesh::Build(tag::grid).shape(tag::segment).start_at(ME).stop_at(MW).divided_in(n_segments(mid_l));
    const Mesh NWM = Mesh::Build(tag::grid).shape(tag::segment).start_at(MW).stop_at(NW).divided_in(n_segments(midW_l));
    const Mesh west = Mesh::Build(tag::grid).shape(tag::segment).start_at(NW).stop_at(SW).divided_in(n_segments(1));
    const Mesh south = Mesh::Build(tag::grid).shape(tag::segment).start_at(SW).stop_at(SE).divided_in(n_segments(1));
    const Mesh east = Mesh::Build(tag::grid).shape(tag::segment).start_at(SE).stop_at(NE).divided_in(n_segments(1));

    const Mesh null_neumann = Mesh::Build(tag::join).meshes({east, west, north_middle}); //podemos ignorar na construção das condições
    const Mesh sources = Mesh::Build(tag::join).meshes({NWM,MNE});
    const Mesh square_boundary = Mesh::Build(tag::join).meshes({south, null_neumann, sources});
    const Mesh north = Mesh::Build(tag::join).meshes({sources, north_middle});

    Mesh inner_boundary = ellipses.total_mesh();
    
    Mesh boundary = Mesh::Build(tag::join).mesh(square_boundary).mesh(inner_boundary);

    const Mesh domain = Mesh::Build(tag::frontal).boundary(boundary).desired_length(h);

    std::map<Cell, size_t> numbering = create_node_numbering(domain, degree);
    
    Eigen::VectorXd solution = build_laplace_solution(base_temp, heat_sources, domain, south, sources, numbering);

    //função objetivo
    double objective = 0;
    double volume = 1 - ellipses.area();

    //Podemos implementar um loop simples mais tarde, integrar trapézios é fácil
    // J(u) = integrate u in north + 10*volume
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

    double final_objective = objective + 10 * volume;

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

SquareSolver::SquareSolver(const double x_m, const double y_m, const double MW_, const double ME_,
    const double heat_sources_, const double base_temp_, const double h,
    const size_t n_ellipses, const bool export_mesh) :
x_max(x_m), y_max(y_m), MW_x(MW_), ME_x(ME_), h(h), num_ellipses(n_ellipses), heat_sources(heat_sources_),
base_temp(base_temp_), export_mesh(export_mesh), ambient(tag::non_existent){
    assert(x_max < 0 || y_max < 0 || //retângulo inválido
        MW_x > ME_x || ME_x > y_max || // MW à direita de ME ou ME à direita de y_max
        num_ellipses < 0);
    
    ambient = Manifold(tag::Euclid, tag::of_dim, 2);
    Function xy = ambient.build_coordinate_system(tag::Lagrange, tag::of_degree, 1);
    Function x = xy[0], y = xy[1];

    //Quadrado Completo (fronteira)
    const Cell NW(tag::vertex); x(NW) = 0; y(NW) = y_max;
    const Cell NE(tag::vertex); x(NE) = x_max; y(NE) = y_max;
    const Cell SW(tag::vertex); x(SW) = 0; y(SW) = 0;
    const Cell SE(tag::vertex); x(SE) = x_max; y(SE) = 0;

    const Cell MW(tag::vertex); x(MW) = MW_x; y(MW) = y_max;
    const Cell ME(tag::vertex); x(ME) = ME_x; y(ME) = y_max;

    double mid_l = ME_x - MW_x;
    auto &midW_l = MW_x;
    double midE_l = 1 - ME_x;
    
    const Mesh MNE = Mesh::Build(tag::grid).shape(tag::segment).start_at(NE).stop_at(ME).divided_in(n_segments(midE_l));
    const Mesh north_middle = Mesh::Build(tag::grid).shape(tag::segment).start_at(ME).stop_at(MW).divided_in(n_segments(mid_l));
    const Mesh NWM = Mesh::Build(tag::grid).shape(tag::segment).start_at(MW).stop_at(NW).divided_in(n_segments(midW_l));
    const Mesh west = Mesh::Build(tag::grid).shape(tag::segment).start_at(NW).stop_at(SW).divided_in(n_segments(1));
    south = Mesh::Build(tag::grid).shape(tag::segment).start_at(SW).stop_at(SE).divided_in(n_segments(1));
    const Mesh east = Mesh::Build(tag::grid).shape(tag::segment).start_at(SE).stop_at(NE).divided_in(n_segments(1));

    const Mesh null_neumann = Mesh::Build(tag::join).meshes({east, west, north_middle}); //podemos ignorar na construção das condições
    sources = Mesh::Build(tag::join).meshes({NWM,MNE});
    square_boundary = Mesh::Build(tag::join).meshes({south, null_neumann, sources});
    north = Mesh::Build(tag::join).meshes({sources, north_middle});
}

double SquareSolver::solve(EllipseBundle& bundle){
    Mesh inner_boundary = bundle.total_mesh();
    
    Mesh boundary = Mesh::Build(tag::join).mesh(square_boundary).mesh(inner_boundary);

    const Mesh domain = Mesh::Build(tag::frontal).boundary(boundary).desired_length(h);

    std::map<Cell, size_t> numbering = create_node_numbering(domain, degree);
    
    Eigen::VectorXd solution = build_laplace_solution(base_temp, heat_sources, domain, south, sources, numbering);

    //função objetivo
    double objective = 0;
    double volume = 1 - bundle.area();

    //Podemos implementar um loop simples mais tarde, integrar trapézios é fácil
    // J(u) = integrate u in north + 10*volume
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

    double final_objective = objective + 10 * volume;

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