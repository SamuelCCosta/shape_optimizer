#include "square_solver.h"
#include <delaunator.hpp>

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


double SquareSolver::solve_frontal(EllipseBundle &bundle) {
    constexpr bool EXPORT_BOUNDARY = true;


    ambient.set_as_working_manifold();
    Function xy = ambient.coordinates();   
    Function x = xy[0], y = xy[1]; 

    Mesh boundary{tag::non_existent};
    
    if (bundle.is_empty()) {
        boundary = square_boundary;
    } else {
        Mesh inner_boundary = bundle.total_mesh();
        boundary = Mesh::Build(tag::join).mesh(square_boundary).mesh(inner_boundary);
    }

    if (EXPORT_BOUNDARY) { boundary.export_to_file(tag::gmsh, "boundary.msh"); }

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

double SquareSolver::solve(EllipseBundle &bundle) {
    constexpr bool EXPORT_BOUNDARY = false;
    constexpr bool DELAUNAY_SOLVER_VERBOSITY = false;

    ambient.set_as_working_manifold();
    Function xy = ambient.coordinates();   
    Function x = xy[0], y = xy[1]; 

    Mesh boundary{tag::non_existent};
    Mesh inner_boundary{tag::non_existent};
    
    if (bundle.is_empty()) {
        boundary = square_boundary;
    } else {
        inner_boundary = bundle.total_mesh();
        boundary = Mesh::Build(tag::join).mesh(square_boundary).mesh(inner_boundary);
    }

    if (EXPORT_BOUNDARY) { boundary.export_to_file(tag::gmsh, "boundary.msh"); }

    // We start with Eigen::Vector2d and we only build the Cell objects when making the final mesh

    // rows of points spread "uniformly" across the entire domain, already cleaned up
    std::vector<Eigen::Vector2d> delaunay_grid = get_delaunay_grid(bundle);

    // get coords and segments vector (coords of boundary + delaunay grid, map with segments)
    std::vector<double> delaunay_coords;
    delaunay_coords.reserve(2 * (delaunay_grid.size() + boundary.number_of(tag::vertices)));
    
    // this way we have a bijection, necessary to build the domain later 
    std::map<Cell, size_t> numbering;
    std::vector<Cell> cells;
    cells.reserve(delaunay_grid.size() + boundary.number_of(tag::vertices));
    size_t counter = 0;
    //empty block to hide iterators
    //we are iterating through the inner boundary first for some filtering later
    {
    if (inner_boundary.exists()) {
        Mesh::Iterator it_inner = inner_boundary.iterator(tag::over_vertices);
        for (it_inner.reset(); it_inner.in_range(); it_inner++) {
            Cell P = *it_inner;
            // could have used create_numbering, but this way we guarantee it has the exact same order
            numbering[P] = counter;
            cells.push_back(P);
            counter++;
            delaunay_coords.push_back(x(P));
            delaunay_coords.push_back(y(P));
        }
    }
    // at this point, counter = n_vertices_inner - 1
    Mesh::Iterator it_outer = square_boundary.iterator(tag::over_vertices);
    for (it_outer.reset(); it_outer.in_range(); it_outer++) {
        Cell P = *it_outer;
        // could have used create_numbering, but this way we guarantee it has the exact same order
        numbering[P] = counter;
        cells.push_back(P);
        counter++;
        delaunay_coords.push_back(x(P));
        delaunay_coords.push_back(y(P));
    }
    // at this point, counter = n_vertices_bdry - 1
    for (auto &point : delaunay_grid) {
        delaunay_coords.push_back(point.x());
        delaunay_coords.push_back(point.y());
        Cell A(tag::vertex, tag::of_coords, {point.x(), point.y()});
        numbering[A] = counter;
        cells.push_back(A);
        counter++;
    }
    // at this point, every single point in the domain has a unique ID
    if (DELAUNAY_SOLVER_VERBOSITY) {
        std::cout << "cells size: " << cells.size() << std::endl;
        std::cout << "numbering size: " <<numbering.size() << std::endl;
        std::cout << "delaunay_coords size: " <<delaunay_coords.size() << std::endl;
    }
    }

    std::map<std::pair<int,int>, Cell> segments;
    auto cell_id = [numbering](Cell A) { return numbering.at(A);};


    auto insert_segment = [cell_id, &segments](Cell seg) {
        Cell A = seg.base().reverse(), B = seg.tip();
        int idA = cell_id(A), idB = cell_id(B);
        int min_id = std::min(idA, idB), max_id = std::max(idA, idB);
        segments.insert({{min_id, max_id}, seg});
    };

    auto get_segment = [cell_id, &segments] (Cell base, Cell tip) {
        int id_base = cell_id(base), id_tip = cell_id(tip);
        int min_id = std::min(id_base, id_tip), max_id = std::max(id_base, id_tip);
        // Check for existence and retrieve in the correct order
        auto finder = segments.find({min_id, max_id});
        if (finder != segments.end()) {
            Cell seg = finder->second;
            return (seg.tip() == tip) ? seg : seg.reverse();
        } else {
            Cell seg(tag::segment, base.reverse(), tip);
            segments.insert({{min_id, max_id}, seg});
            return seg;
        }
    };

    // Filling up segments map
    {
    Mesh::Iterator it = boundary.iterator(tag::over_segments);
    for (it.reset(); it.in_range(); it++) { insert_segment(*it); }
    }
    // Delaunay triangulation
    delaunator::Delaunator delaunay(delaunay_coords);
    std::vector<std::size_t> triangles = delaunay.triangles;
    
    if (DELAUNAY_SOLVER_VERBOSITY) {
        std::cout << "number of triangles before filtering: " << triangles.size() / 3 << std::endl;
    }

    // build final mesh
    Mesh domain(tag::fuzzy, tag::of_dim, 2);
    size_t bdry_upper_bound = inner_boundary.exists() ? inner_boundary.number_of(tag::vertices) : 0; //all ellipse boundary points have id < this
    for (size_t i = 0; i < triangles.size(); i += 3) {
        // Delaunator outputs triangles in clockwise order for a Y-up coordinate system.
        // maniFEM expects counter-clockwise order for positive triangles, so we swap id1 and id2.
        size_t id0 = triangles[i], id1 = triangles[i+2], id2 = triangles[i+1];
        // Check if the triangle has every vertice in the ellipse boundary; if so, skip the current iteration
        if (id0 < bdry_upper_bound && id1 < bdry_upper_bound && id2 < bdry_upper_bound) { continue; }

        // WARNING: IT MIGHT REMOVE TRIANGLES WITH POINTS THAT LIE IN DIFFERENT ELLIPSES (don't have a fix, for now)
        Cell triangle(tag::triangle, get_segment(cells[id0], cells[id1]), get_segment(cells[id1], cells[id2]), get_segment(cells[id2], cells[id0]));
        
        triangle.add_to(domain);
    }

    if (DELAUNAY_SOLVER_VERBOSITY) {
        std::cout << "number of triangles after filtering: " << domain.number_of(tag::cells_of_max_dim) << std::endl;
    }
    
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


std::vector<Eigen::Vector2d> SquareSolver::get_delaunay_grid(EllipseBundle &bundle) {
    std::vector<Eigen::Vector2d> points;
    // delta calculation to reject points too near to the ellipses
    std::vector<double> deltas;
    deltas.reserve(bundle.bundle.size());

    const double threshold = 0.4 * h;
    for (const auto &ellipse : bundle.bundle) {
        double invs_ellipse_minor = 1 / ellipse.parametrize_matrix.col(1).norm();
        double delta = 2 * threshold * invs_ellipse_minor + threshold * threshold * invs_ellipse_minor * invs_ellipse_minor;
        deltas.push_back(delta);
    }

    const double sq3_over_2 = std::sqrt(3.0)/2;
    int n_vertical_divisions =  static_cast<int>(std::round(x_max / h * sq3_over_2));
    double dy = y_max / n_vertical_divisions;
    // i = 1, ..., n_vert_divs - 1
    for (int i = n_vertical_divisions - 1; i > 0; i--) {
        double y = i * dy;
        double offset = (i % 2) * h/2;
        int j_max = static_cast<int>(std::round( (x_max - offset) / h ));

        for (int j = 1; j < j_max - 1; j++){
            double x = offset + j * h;
            Eigen::Vector2d point(x, y);
            bool is_outside_all = true;
            int k = 0;
            for (const auto &ellipse : bundle.bundle) {
                double delta = deltas[k];
                // in the worst case, a point can be excluded at approx 2.72*threshold (depends on eccentricity, with the formula 1/sqrt(1-ecc^2))
                if (ellipse.evaluate_at(point) < 1.0 + delta) { //slightly bigger to exclude points close to boundary
                    is_outside_all = false;
                    break;
                }
                k++;
            }
            if (is_outside_all) {
                points.push_back(point);
            }
        }
    }

    return points;
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
	for ( size_t j = 0; j < size_matrix; j++ ) {
		if ( matrix_A .coeff ( cell_id, j ) != 0. ) {
			matrix_A .coeffRef ( cell_id, j ) = 0.;
        }
    }

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
    
    Mesh::Iterator it = sources.iterator(tag::over_cells_of_max_dim);
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