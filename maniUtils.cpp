#include "maniUtils.h"
#include <iostream>

using namespace maniFEM;


//No caso de condições de fronteira de Dirichlet
void impose_value_of_unknown
( Eigen::SparseMatrix < double >& matrix_A, Eigen::VectorXd& vector_b,
  const size_t i, const double val)
{	assert ( ! matrix_A .IsRowMajor );
	size_t size_matrix = matrix_A .innerSize();

	// apagar linha i
	for ( size_t j = 0; j < size_matrix; j++ )
		if ( matrix_A .coeff ( i, j ) != 0. )
			matrix_A .coeffRef ( i, j ) = 0.;

	// apagar coluna i e mudar vetor
	for ( size_t j = 0; j < size_matrix; j++ )
	{	if ( matrix_A .coeff ( j, i ) == 0. )  continue;
		double & Aji = matrix_A .coeffRef ( j, i );
		vector_b (j) -= Aji * val;
		Aji = 0.;                                         }

	// colocar a_ii = 1 e b_i = val
	matrix_A .coeffRef ( i, i ) = 1.;
	vector_b (i) = val;                                     }


std::pair<Mesh, Mesh> create_mesh_from_implicit(Manifold& ambient_manifold, 
    const Function& implicit_eq, const double h, const Cell& start_point){

        Manifold mani_implicit = ambient_manifold.implicit(implicit_eq == 0); //conjunto de nível 0

        Mesh bdry = Mesh::Build(tag::frontal).entire_manifold().start_at(start_point).desired_length(h);

        ambient_manifold.set_as_working_manifold();
        Mesh final_mesh = Mesh::Build (tag::frontal). boundary(bdry). desired_length(h);

        return{final_mesh, bdry};
}

std::map<Cell, size_t> create_node_numbering(const Mesh& mesh, const size_t degree = 1){
    std::map<Cell, size_t> numbering;
    size_t counter = 0;

    if (degree == 1) {
    Mesh::Iterator it = mesh.iterator(tag::over_vertices);
    for (it .reset() ; it .in_range(); it++){
        Cell V = *it; 
        numbering [V] = counter; 
        counter++;}
    assert(counter == numbering.size());
    } else if (degree == 2) {
    Mesh::Iterator it = mesh.iterator(tag::over_segments);
    for (it.reset() ; it.in_range(); it++){
        Cell seg = *it; 
        numbering [seg] = counter; 
        counter++;}
    assert(counter == numbering.size());
    } else {
        std::cerr << "Not implemented, degree must be 1 or 2" << std::endl;
    };

    return numbering;
}

