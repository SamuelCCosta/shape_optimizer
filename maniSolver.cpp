#include "maniSolver.h"
#include "maniUtils.h"

Eigen::VectorXd build_poisson_solution(const Function& f, const Function& boundary_value, const Function& normal_deriv,
    const Mesh& domain, const Mesh& dirichlet_boundary, const Mesh& neumann_boundary, const Mesh& boundary,
    const std::map<Cell, size_t>& numbering, const size_t& degree){
    
    assert(degree == 1 || degree == 2);
    
    Manifold& active_manifold = Manifold::working;
    Function xy = active_manifold.coordinates();
    Function x = xy[0], y = xy[1];

    std::optional<FiniteElement> fe_opt;
    std::optional<Integrator> integr_opt;
    std::optional<FiniteElement> fe_bdry_opt;
    std::optional<Integrator> integr_bdry_opt;

    if (degree == 1) {
        fe_opt.emplace(tag::with_master, tag::triangle, tag::Lagrange, tag::of_degree, 1);
        integr_opt.emplace(fe_opt -> set_integrator(tag::Gauss, tag::tri_6));
    } else if (degree == 2) {
        fe_opt.emplace(tag::with_master, tag::triangle, tag::Lagrange, tag::of_degree, 2, tag::straight);
        integr_opt.emplace(fe_opt -> set_integrator(tag::Gauss, tag::tri_6));    
    }
    // Fronteira é sempre grau 1
    FiniteElement fe_bdry(tag::with_master, tag::segment, tag::Lagrange, tag::of_degree, 1);
    Integrator integr_bdry = fe_bdry.set_integrator(tag::Gauss, tag::seg_3);
    
    FiniteElement& fe = fe_opt.value();

    //Definir matriz A e vetor b
    size_t size_matrix = numbering.size();
    Eigen::SparseMatrix<double> matrix_A(size_matrix, size_matrix);
    Eigen::VectorXd vector_b(size_matrix); vector_b.setZero();

    //reservar espaço dependendo da complexidade
    if (degree == 1) {
        matrix_A.reserve(Eigen::VectorXi::Constant(size_matrix, 8));
    } else if (degree == 2) {
        matrix_A.reserve(Eigen::VectorXi::Constant(size_matrix, 16));
    }

    if (degree == 1) {
        Mesh::Iterator it = domain.iterator(tag::over_cells_of_max_dim);
        for(it.reset(); it.in_range(); it++){
            Cell tri = *it;
            fe.dock_on(tri); 
            //Iterar duas vezes sob vértices
            Mesh::Iterator itV = tri.boundary().iterator(tag::over_vertices);
            Mesh::Iterator itW = tri.boundary().iterator(tag::over_vertices);

            for(itV.reset(); itV.in_range(); itV++){
                Cell V = *itV;
                Function phi_V = fe.basis_function(V),
                            phi_V_dx = phi_V.deriv(x),
                            phi_V_dy = phi_V.deriv(y);
            
                for(itW.reset(); itW.in_range(); itW++){
                    Cell W = *itW;

                    Function phi_W = fe.basis_function(W),
                                phi_W_dx = phi_W.deriv(x),
                                phi_W_dy = phi_W.deriv(y);
                    
                    //Integrar a(u,v) e adicionar à matriz
                    matrix_A.coeffRef(numbering.at(V), numbering.at(W)) +=
                        fe.integrate(phi_V_dx * phi_W_dx + phi_V_dy * phi_W_dy);
                }

                vector_b(numbering.at(V)) += fe.integrate(f * phi_V);
            }
        } 
    } else if (degree == 2) {
        Mesh::Iterator it = domain.iterator(tag::over_cells_of_max_dim);
        for(it.reset(); it.in_range(); it++){
            Cell tri = *it;
            fe.dock_on(tri);
            
            Mesh::Iterator it_V = tri.boundary().iterator(tag::over_vertices);
            Mesh::Iterator it_W = tri.boundary().iterator(tag::over_vertices);
            for (it_V.reset(); it_V.in_range(); it_V++){
                Cell V = *it_V;

                Function psi_V = fe .basis_function(V),
                             psi_V_dx = psi_V.deriv(x),
                             psi_V_dy = psi_V.deriv(y);
                Cell seg_V = tri.boundary().cell_in_front_of(V);
                Function psi_sV = fe.basis_function(seg_V),
                               psi_sV_dx = psi_sV.deriv(x),
                               psi_sV_dy = psi_sV.deriv(y);
                Cell sV = seg_V.get_positive();

                for (it_W.reset(); it_W.in_range(); it_W++){
                    Cell W = *it_W;

                    Function psi_W = fe.basis_function(W),
                                psi_W_dx = psi_W.deriv(x),
                                psi_W_dy = psi_W.deriv(y);
                    Cell seg_W = tri.boundary().cell_in_front_of(W);
                    Function psi_sW = fe.basis_function(seg_W),
                                   psi_sW_dx = psi_sW.deriv(x),
                                   psi_sW_dy = psi_sW.deriv(y);
                    Cell sW = seg_W.get_positive();

                    matrix_A.coeffRef (numbering.at(V), numbering.at(W) ) +=
                        fe.integrate (psi_V_dx * psi_W_dx + psi_V_dy * psi_W_dy);
                    matrix_A.coeffRef (numbering.at(V), numbering.at(sW) ) +=
                        fe.integrate (psi_V_dx * psi_sW_dx + psi_V_dy * psi_sW_dy);
                    matrix_A.coeffRef (numbering.at(sV), numbering.at(W) ) +=
                        fe.integrate (psi_sV_dx * psi_W_dx + psi_sV_dy * psi_W_dy);
                    matrix_A.coeffRef ( numbering.at(sV), numbering.at(sW) ) +=
                        fe.integrate (psi_sV_dx * psi_sW_dx + psi_sV_dy * psi_sW_dy);
                }

                vector_b(numbering.at(V)) += fe.integrate(f * psi_V);
                vector_b(numbering.at(sV)) += fe.integrate(f * psi_sV);
            }
        }
    }

    // Condição de Neumann
    if (degree == 1) {
        Mesh::Iterator it = neumann_boundary.iterator(tag::over_cells_of_max_dim);
        for(it.reset(); it.in_range(); it++){
            Cell seg = *it;
            fe_bdry.dock_on(seg);
            Mesh::Iterator it_vert = seg.boundary().iterator(tag::over_vertices, tag::force_positive);
            for(it_vert.reset(); it_vert.in_range(); it_vert++){
                Cell V = *it_vert;
                Function phi_V = fe_bdry.basis_function(V);
                vector_b(numbering.at(V)) +=
                    fe_bdry.integrate(normal_deriv * phi_V);
            }
        }
    } else if (degree == 2) {
        Mesh::Iterator it = neumann_boundary.iterator(tag::over_cells_of_max_dim);
        for(it.reset(); it.in_range(); it++){
            Cell seg = *it;
            fe_bdry.dock_on(seg);
            Mesh::Iterator it_vert = seg.boundary().iterator(tag::over_vertices, tag::force_positive);
            for(it_vert.reset(); it_vert.in_range(); it_vert++){
                Cell V = *it_vert;
                Function phi_V = fe_bdry.basis_function(V);
                vector_b(numbering.at(V)) +=
                    fe_bdry.integrate(normal_deriv * phi_V);
            }
            // Um dia haverão elementos finitos de grau 2 para segmentos
            /*Function phi_seg = fe_bdry.basis_function(seg);
            Cell s = seg.get_positive();
            vector_b(numbering.at(s)) +=
                fe_bdry.integrate(normal_deriv * phi_seg);*/
        }
    }

    //Condição de Dirichlet
    if (degree == 1) {
        Mesh::Iterator it = dirichlet_boundary. iterator(tag::over_vertices);
        for(it.reset(); it.in_range(); it++){
            Cell P = *it;
            size_t i = numbering.at(P);
            //Forçar valor dos nodos da fronteira para boundary_value
            impose_value_of_unknown(matrix_A, vector_b, i, boundary_value(P));
        }
    } else if (degree == 2){
        Mesh::Iterator it = dirichlet_boundary. iterator(tag::over_vertices);
        for(it.reset(); it.in_range(); it++){
            Cell P = *it;
            size_t i = numbering.at(P);
            impose_value_of_unknown(matrix_A, vector_b, i, boundary_value(P));

            Cell seg = boundary.cell_in_front_of(P, tag::surely_exists);
            Cell Q = seg .tip();
            i = numbering.at(seg.get_positive());
            Cell midpoint(tag::vertex); x(midpoint) = (x(P)+x(Q))/2.; y(midpoint) = (y(P)+y(Q))/2.;
            impose_value_of_unknown(matrix_A, vector_b, i, boundary_value(midpoint));
        }
    }


    //Resolver sistema linear
    Eigen::ConjugateGradient <Eigen::SparseMatrix<double>,
                            Eigen::Lower | Eigen::Upper> cg;

    cg.compute(matrix_A);

    Eigen::VectorXd solution = cg.solve(vector_b);
    if(cg.info() != Eigen::Success)
        std::cout << "Eigen solver failed" << std::endl;

    return solution;
    }


Eigen::VectorXd build_poisson_solution(const Function& f, const Function& boundary_value, const Function& normal_deriv,
    const Mesh& domain, const Mesh& dirichlet_boundary, const Mesh& neumann_boundary, const std::map<Cell, size_t>& numbering, const tag::HandCoded&){

    Manifold& active_manifold = Manifold::working;
    Function xy = active_manifold.coordinates();
    Function x = xy[0], y = xy[1];

    FiniteElement fe(tag::with_master, tag::triangle, tag::Lagrange, tag::of_degree, 1);
    Integrator integr = fe.set_integrator(tag::Gauss, tag::tri_6);

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
        fe_hand.dock_on(tri); //problema aqui
        fe.dock_on(tri);
        //Iterar duas vezes sob vértices
        Mesh::Iterator itV = tri.boundary().iterator(tag::over_vertices);
        Mesh::Iterator itW = tri.boundary().iterator(tag::over_vertices);
    
        for(itV.reset(); itV.in_range(); itV++){
            Cell V = *itV;
            Function phi_V = fe_hand.basis_function(V);
            Function nothand_phi_V = fe.basis_function(V);

            for(itW.reset(); itW.in_range(); itW++){
                Cell W = *itW;
                Function phi_W = fe_hand.basis_function(W);

                std::vector<double> result = fe_hand.integrate(tag::pre_computed,
                                              tag::replace, bf1, tag::by, phi_V,
                                              tag::replace, bf2, tag::by, phi_W);

                matrix_A.coeffRef(numbering.at(V),numbering.at(W)) += result[0];
            }
            
            vector_b(numbering.at(V)) += fe.integrate(f * nothand_phi_V);
        }
    }
    }   
    //Condição de Neumann
    {
    Mesh::Iterator it = neumann_boundary.iterator(tag::over_cells_of_max_dim);
    for(it.reset(); it.in_range(); it++){
        Cell seg = *it;
        fe_bdry.dock_on(seg);
        Mesh::Iterator it_vert = seg.boundary().iterator(tag::over_vertices, tag::force_positive);
        for(it_vert.reset(); it_vert.in_range(); it_vert++){
            Cell V = *it_vert;
            Function phi_V = fe_bdry.basis_function(V);
            vector_b(numbering.at(V)) +=
                fe_bdry.integrate(normal_deriv * phi_V);
        }
    }
    }
    //Condição de Dirichlet
    {
    Mesh::Iterator it = dirichlet_boundary.iterator(tag::over_vertices);
    for(it.reset(); it.in_range(); it++){
        Cell P = *it;
        size_t i = numbering.at(P);
        //Forçar valor dos nodos da fronteira para boundary_value
        impose_value_of_unknown(matrix_A, vector_b, i, boundary_value(P));
    }
    }
    //Resolver sistema linear
    Eigen::ConjugateGradient <Eigen::SparseMatrix<double>,
                            Eigen::Lower | Eigen::Upper> cg;

    cg.compute(matrix_A);

    Eigen::VectorXd solution = cg.solve(vector_b);
    if(cg.info() != Eigen::Success)
        std::cout << "Eigen solver failed" << std::endl;

    return solution;
}

//hand coded Lagrange P1 finite elements
Eigen::VectorXd build_laplace_solution(const Function& boundary_value, const Function& normal_deriv, const Mesh& domain,
    const Mesh& dirichlet_boundary, const Mesh& neumann_boundary, const std::map<Cell, size_t>& numbering){
        
    Manifold& active_manifold = Manifold::working;
    Function xy = active_manifold.coordinates();
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
        Mesh::Iterator itW = tri.boundary().iterator(tag::over_vertices);
    
        for(itV.reset(); itV.in_range(); itV++){
            Cell V = *itV;
            Function phi_V = fe_hand.basis_function(V);

            for(itW.reset(); itW.in_range(); itW++){
                Cell W = *itW;
                Function phi_W = fe_hand.basis_function(W);

                std::vector<double> result = fe_hand.integrate(tag::pre_computed,
                                              tag::replace, bf1, tag::by, phi_V,
                                              tag::replace, bf2, tag::by, phi_W);

                matrix_A.coeffRef(numbering.at(V),numbering.at(W)) += result[0];
            }
        }
    }
    }   
    //Condição de Neumann
    {
    Mesh::Iterator it = neumann_boundary.iterator(tag::over_cells_of_max_dim);
    for(it.reset(); it.in_range(); it++){
        Cell seg = *it;
        fe_bdry.dock_on(seg);
        Mesh::Iterator it_vert = seg.boundary().iterator(tag::over_vertices, tag::force_positive);
        for(it_vert.reset(); it_vert.in_range(); it_vert++){
            Cell V = *it_vert;
            Function phi_V = fe_bdry.basis_function(V);
            vector_b(numbering.at(V)) +=
                fe_bdry.integrate(normal_deriv * phi_V);
        }
    }
    }
    //Condição de Dirichlet
    {
    Mesh::Iterator it = dirichlet_boundary.iterator(tag::over_vertices);
    for(it.reset(); it.in_range(); it++){
        Cell P = *it;
        size_t i = numbering.at(P);
        //Forçar valor dos nodos da fronteira para boundary_value
        impose_value_of_unknown(matrix_A, vector_b, i, boundary_value(P));
    }
    }
    //Resolver sistema linear
    Eigen::ConjugateGradient <Eigen::SparseMatrix<double>,
                            Eigen::Lower | Eigen::Upper> cg;

    cg.compute(matrix_A);

    Eigen::VectorXd solution = cg.solve(vector_b);
    if(cg.info() != Eigen::Success)
        std::cout << "Eigen solver failed" << std::endl;

    return solution;
}

//para condições de fronteira constantes
Eigen::VectorXd build_laplace_solution(const double boundary_value, const double normal_deriv, const Mesh& domain,
    const Mesh& dirichlet_boundary, const Mesh& neumann_boundary, const std::map<Cell, size_t>& numbering){
        
    Manifold& active_manifold = Manifold::working;
    Function xy = active_manifold.coordinates();
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
        fe_hand.dock_on(tri); //problema aqui
        //Iterar duas vezes sob vértices
        Mesh::Iterator itV = tri.boundary().iterator(tag::over_vertices);
        Mesh::Iterator itW = tri.boundary().iterator(tag::over_vertices);
    
        for(itV.reset(); itV.in_range(); itV++){
            Cell V = *itV;
            Function phi_V = fe_hand.basis_function(V);

            for(itW.reset(); itW.in_range(); itW++){
                Cell W = *itW;
                Function phi_W = fe_hand.basis_function(W);

                std::vector<double> result = fe_hand.integrate(tag::pre_computed,
                                              tag::replace, bf1, tag::by, phi_V,
                                              tag::replace, bf2, tag::by, phi_W);

                matrix_A.coeffRef(numbering.at(V),numbering.at(W)) += result[0];
            }
        }
    }
    }   
    //Condição de Neumann
    {
    Mesh::Iterator it = neumann_boundary.iterator(tag::over_cells_of_max_dim);
    for(it.reset(); it.in_range(); it++){
        Cell seg = *it;
        fe_bdry.dock_on(seg);
        Mesh::Iterator it_vert = seg.boundary().iterator(tag::over_vertices, tag::force_positive);
        for(it_vert.reset(); it_vert.in_range(); it_vert++){
            Cell V = *it_vert;
            Function phi_V = fe_bdry.basis_function(V);
            vector_b(numbering.at(V)) +=
                normal_deriv * fe_bdry.integrate(phi_V);
        }
    }
    }
    //Condição de Dirichlet
    {
    Mesh::Iterator it = dirichlet_boundary.iterator(tag::over_vertices);
    for(it.reset(); it.in_range(); it++){
        Cell P = *it;
        size_t i = numbering.at(P);
        //Forçar valor dos nodos da fronteira para boundary_value
        impose_value_of_unknown(matrix_A, vector_b, i, boundary_value);
    }
    }
    //Resolver sistema linear
    Eigen::ConjugateGradient <Eigen::SparseMatrix<double>,
                            Eigen::Lower | Eigen::Upper> cg;

    cg.compute(matrix_A);

    Eigen::VectorXd solution = cg.solve(vector_b);
    if(cg.info() != Eigen::Success)
        std::cout << "Eigen solver failed" << std::endl;

    return solution;
}