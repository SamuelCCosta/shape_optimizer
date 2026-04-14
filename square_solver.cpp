#include "square_solver.h"

SquareSolver::SquareSolver(std::map<std::string, double> &geometric_config,
    const double h_param,
    const double heat_srcs,
    const double base_tmp,
    const bool export_dom,
    const bool export_res) :
    x_max(geometric_config["x_max"]), y_max(geometric_config["y_max"]), MW_x(geometric_config["MW_x"]),
    ME_x(geometric_config["ME_x"]), h(h_param), heat_sources(heat_srcs), base_temp(base_tmp),
    export_domain(export_dom), export_result(export_res)
    {
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
    double &midW_l = MW_x;
    double midE_l = x_max - ME_x;

    const Mesh MNE = Mesh::Build(tag::grid).shape(tag::segment).start_at(NE).stop_at(ME)
                                                            .divided_in(n_segments(midE_l));
    const Mesh north_middle = Mesh::Build(tag::grid).shape(tag::segment).start_at(ME).stop_at(MW)
                                                            .divided_in(n_segments(mid_l));
    const Mesh NWM = Mesh::Build(tag::grid).shape(tag::segment).start_at(MW).stop_at(NW)
                                                            .divided_in(n_segments(midW_l));
    const Mesh west = Mesh::Build(tag::grid).shape(tag::segment).start_at(NW).stop_at(SW)
                                                            .divided_in(n_segments(y_max));
    south = Mesh::Build(tag::grid).shape(tag::segment).start_at(SW).stop_at(SE)
                                                            .divided_in(n_segments(x_max));
    const Mesh east = Mesh::Build(tag::grid).shape(tag::segment).start_at(SE).stop_at(NE)
                                                            .divided_in(n_segments(y_max));
    
    //podemos ignorar na construção das condições
    const Mesh null_neumann = Mesh::Build(tag::join).meshes({east, west, north_middle});

    sources = Mesh::Build(tag::join).meshes({NWM,MNE});
    square_boundary = Mesh::Build(tag::join).meshes({south, null_neumann, sources});
    north = Mesh::Build(tag::join).meshes({sources, north_middle});
}

SquareSolver::SquareSolver(std::map<std::string, double> &geometric_config,
    const double h_param,
    const double heat_srcs,
    const double base_tmp) :
    SquareSolver(geometric_config, h_param, heat_srcs, base_tmp, false, false) {}


double SquareSolver::solve(EllipseBundle &bundle) {
    ambient.set_as_working_manifold();
    Function xy = ambient.coordinates();   
    Function x = xy[0], y = xy[1]; 

    Mesh inner_boundary = bundle.total_mesh();
    Mesh boundary = Mesh::Build(tag::join).mesh(square_boundary).mesh(inner_boundary);
    boundary.export_to_file(tag::gmsh, "boundary_debug.msh"); //debug

    const Mesh domain = Mesh::Build(tag::frontal).boundary(boundary).desired_length(h);
    std::map<Cell, size_t> numbering = create_numbering(domain);

    if (export_domain) { domain.export_to_file(tag::gmsh, "domain.msh"); }

    Eigen::VectorXd solution = build_laplace_solution(domain, numbering);

    if (export_result) {
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

    return objective_no_penalty(solution, numbering); 
}

std::map<Cell, size_t> SquareSolver::create_numbering(const Mesh& mesh) {
    std::map<Cell, size_t> numbering;
    size_t counter = 0;

    Mesh::Iterator it = mesh.iterator(tag::over_vertices);
    for (it .reset() ; it .in_range(); it++){
        Cell V = *it; 
        numbering [V] = counter; 
        counter++;}
    assert(counter == numbering.size());

    return numbering;
}

Eigen::VectorXd SquareSolver::build_laplace_solution(const Mesh &domain, const std::map<Cell,size_t> &numbering) {
    Function xy = ambient.coordinates();
    Function x = xy[0], y = xy[1];

    FiniteElement fe_hand(tag::triangle, tag::Lagrange, tag::of_degree, 1);
    Integrator hand_integr = fe_hand.set_integrator(tag::hand_coded);
    Function bf1(tag::basis_function, tag::within, fe_hand),
             bf2(tag::basis_function, tag::within, fe_hand);
    
    fe_hand.pre_compute(tag::for_given, tag::basis_functions, bf1, bf2, tag::integral_of,
                            {bf1 .deriv(x) * bf2 .deriv(x) + bf1 .deriv(y) * bf2 .deriv(y)});

    FiniteElement fe_bdry(tag::with_master, tag::segment, tag::Lagrange, tag::of_degree, 1);
    Integrator integr_bdry= fe_bdry.set_integrator(tag::Gauss, tag::seg_3);
    
    //Definir matriz A e vetor b
    size_t size_matrix = numbering.size();
    Eigen::SparseMatrix<double> matrix_A(size_matrix, size_matrix);
    Eigen::VectorXd vector_b(size_matrix); vector_b.setZero();

    matrix_A.reserve(Eigen::VectorXi::Constant(size_matrix, 8));


    //Construção do sistema
    {
    Mesh::Iterator it = domain.iterator(tag::over_cells_of_max_dim);
    for(it.reset(); it.in_range(); it++){
        Cell tri = *it;
        fe_hand.dock_on(tri);
        //Iterar duas vezes sob vértices
        Mesh::Iterator itV = tri.boundary().iterator(tag::over_vertices);

        for(itV.reset(); itV.in_range(); itV++){
            Cell V = *itV;
            Function phi_V = fe_hand.basis_function(V);
            Mesh::Iterator itW = itV;

            for(; itW.in_range(); itW++){
                Cell W = *itW;
                Function phi_W = fe_hand.basis_function(W);

                std::vector<double> result = fe_hand.integrate(tag::pre_computed,
                                              tag::replace, bf1, tag::by, phi_V,
                                              tag::replace, bf2, tag::by, phi_W);

                size_t idxV = numbering.at(V), idxW = numbering.at(W); 
                matrix_A.coeffRef(idxV, idxW) += result[0];
                //Check if it isn't in the diagonal
                if (idxW != idxV) {matrix_A.coeffRef(idxW, idxV) += result[0];}
            }
        }
    }
    }

    //Condição de Neumann
    {
    Mesh::Iterator it = sources.iterator(tag::over_cells_of_max_dim);
    for(it.reset(); it.in_range(); it++){
        Cell seg = *it;
        fe_bdry.dock_on(seg);
        Mesh::Iterator it_vert = seg.boundary().iterator(tag::over_vertices, tag::force_positive);
        for(it_vert.reset(); it_vert.in_range(); it_vert++){
            Cell V = *it_vert;
            Function phi_V = fe_bdry.basis_function(V);
            vector_b(numbering.at(V)) +=
                heat_sources * fe_bdry.integrate(phi_V);
        }
    }
    }

    //Condição de Dirichlet
    {
    Mesh::Iterator it = south.iterator(tag::over_vertices);
    for(it.reset(); it.in_range(); it++){
        Cell P = *it;
        size_t i = numbering.at(P);
        //Forçar valor dos nodos da fronteira para base_temp
        impose_value_of_unknown(matrix_A, vector_b, i, base_temp);
    }
    }

    //Resolver sistema linear
    Eigen::ConjugateGradient <Eigen::SparseMatrix<double>,
                            Eigen::Lower | Eigen::Upper> cg;

    cg.compute(matrix_A);

    Eigen::VectorXd solution = cg.solve(vector_b);
    if(cg.info() != Eigen::Success) {
        std::cout << "Eigen solver failed" << std::endl;
    }

    return solution;
}

void SquareSolver::impose_value_of_unknown(Eigen::SparseMatrix<double> &matrix_A, 
    Eigen::VectorXd &vector_b, const size_t cell_id, const double val) {
    assert ( ! matrix_A .IsRowMajor );
	size_t size_matrix = matrix_A .innerSize();

	// apagar linha i (cell_id = i)
	for ( size_t j = 0; j < size_matrix; j++ )
		if ( matrix_A .coeff ( cell_id, j ) != 0. )
			matrix_A .coeffRef ( cell_id, j ) = 0.;

	// apagar coluna i e mudar vetor
	for ( size_t j = 0; j < size_matrix; j++ )
	{	if ( matrix_A .coeff ( j, cell_id ) == 0. )  continue;
		double & Aji = matrix_A .coeffRef ( j, cell_id );
		vector_b (j) -= Aji * val;
		Aji = 0.;                                         }

	// colocar a_ii = 1 e b_i = val
	matrix_A .coeffRef ( cell_id, cell_id ) = 1.;
	vector_b (cell_id) = val;
}

double SquareSolver::objective_no_penalty(const Eigen::VectorXd &solution, const std::map<Cell,size_t> &numbering) {
    double objective = 0;

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
    return objective;
}