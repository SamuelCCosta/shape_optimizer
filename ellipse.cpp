#include "ellipse.h"

Ellipse::Ellipse(double x, double y, double A, double B, double C) {
    this-> quadratic_form << A, B, B, C;
    this-> center << x, y;
    double det = quadratic_form.determinant();
    if (A < 0 || det < 0) {
        throw std::invalid_argument("The Matrix is not positive definite");
    }

    this-> bounds << std::sqrt(C/det), std::sqrt(A/det); //width, height 

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> solver(quadratic_form);

    Eigen::Vector2d eigenvalues = solver.eigenvalues(); //sorted in increasing order
    Eigen::Matrix2d eigenvector_matrix = solver.eigenvectors();

    Eigen::Matrix2d scaling = Eigen::Matrix2d::Zero();
    scaling(0,0) = 1.0 / std::sqrt(eigenvalues(0));
    scaling(1,1) = 1.0 / std::sqrt(eigenvalues(1));

    this-> parametrize_matrix = eigenvector_matrix * scaling;

    double eccentricity_2 = 1.0 - eigenvalues(0)/eigenvalues(1);
    if (eccentricity_2 > 0.93 * 0.93) { throw std::invalid_argument("Eccentricity > 0.93"); }
}

Mesh Ellipse::get_mesh(const double h) const {
    Manifold RR2 = Manifold::working;
    Function xy = RR2.coordinates();
    Function x = xy[0], y = xy[1];

    double xc = center[0], yc = center[1];
    double A = quadratic_form(0,0), B = quadratic_form(1,0), C = quadratic_form(1,1);
    Function implicit_eq = A * (x-xc) * (x-xc) + 2 * B * (x-xc) * (y-yc) + C * (y-yc) * (y-yc); 

    Manifold representation = RR2.implicit(implicit_eq == 1);

    //starting point and starting derivative
    Eigen::Vector2d start_point = point_at(0), start_derivative = - derivative_at(0);
    Cell start(tag::vertex, tag::of_coordinates, {start_point.x(), start_point.y()});
    std::vector<double> direction = {start_derivative.x(), start_derivative.y()};

    Mesh mesh = Mesh::Build(tag::frontal).entire_manifold().start_at(start).towards(direction).desired_length(h);

    RR2.set_as_working_manifold();

    return mesh;
}

bool EllipseBundle::intersects(const Ellipse &e1, const Ellipse &e2) const {
    // Inexpensive: check if centers are inside
    if (e2.evaluate_at(e1.center) <= 1.0) return true;
    if (e1.evaluate_at(e2.center) <= 1.0) return true;

    // Semi-Inexpensive: bounding boxes
    Eigen::Vector2d h_vec(h,h);

    Eigen::Vector2d e1_bottom_left = e1.center - e1.bounds;
    Eigen::Vector2d e1_top_right = e1.center + e1.bounds;

    Eigen::Vector2d e2_bottom_left = e2.center - e2.bounds;
    Eigen::Vector2d e2_top_right = e2.center + e2.bounds;

    //vec1 > vec2 iff at least one of the components is bigger
    if (is_greater(e2_bottom_left, e1_top_right + h_vec) || 
        is_greater(e1_bottom_left, e2_top_right + h_vec))
        { return false; }

    return robust_intersect(e1,e2);
}

bool EllipseBundle::robust_intersect(const Ellipse &e1, const Ellipse &e2) const {
    const double h_sq = h * h;
    const double eps = 1e-10;
    const int max_iterations = 10;
    
    auto [theta1, theta2] = starting_parameters(e1, e2);
    
    for (int i = 0; i < max_iterations; i++) {
        // Useful numbers
        double cos_1 = std::cos(theta1), sin_1 = std::sin(theta1);
        double cos_2 = std::cos(theta2), sin_2 = std::sin(theta2);

        Eigen::Vector2d P1 = e1.point_at(Eigen::Vector2d(cos_1, sin_1));
        Eigen::Vector2d P2 = e2.point_at(Eigen::Vector2d(cos_2, sin_2));

        // Sanity check: are the points outside of the other ellipse?
        if (e2.is_inside(P1) || e1.is_inside(P2)) return true;

        // Check if the distance^2 between points is less than h^2
        Eigen::Vector2d difference = P1 - P2;
        if (difference.squaredNorm() <= h_sq) return true;

        Eigen::Vector2d P_prime_1 = e1.derivative_at(Eigen::Vector2d(- sin_1, cos_1));
        Eigen::Vector2d P_prime_2 = e2.derivative_at(Eigen::Vector2d(- sin_2, cos_2));

        // (dot_1,dot_2) are the functions we are trying to equal to zero
        double dot_1 = difference.dot(P_prime_1), dot_2 = difference.dot(P_prime_2);
        if (std::abs(dot_1) < eps && std::abs(dot_2) < eps) break; //optimization is complete

        // Update theta1 and theta2
        Eigen::Vector2d P_doubleprime_1 = e1.second_derivative_at(Eigen::Vector2d(- cos_1, - sin_1));
        Eigen::Vector2d P_doubleprime_2 = e2.second_derivative_at(Eigen::Vector2d(- cos_2, - sin_2));

        // Hessian Matrix
        double H00 = P_prime_1.squaredNorm() + difference.dot(P_doubleprime_1);
        double H01 = - P_prime_1.dot(P_prime_2);
        double H11 = P_prime_2.squaredNorm() - difference.dot(P_doubleprime_2);

        double hessian_determinant = H00 * H11 - H01 * H01;
        if (std::abs(hessian_determinant) < eps) return true; //Degenerate Hessian

        // Solve linear system with Cramer's Rule
        // [H00 H01] [ delta_1 ] = [ - dot_1 ]
        // [H01 H11] [ delta_2 ]   [ - dot_2 ]

        double delta_1 = (- dot_1 * H11 + H01 * dot_2) / hessian_determinant;
        double delta_2 = (H00 * (- dot_2) + dot_1 * H01) / hessian_determinant;

        // update parameters
        theta1 += delta_1, theta2 += delta_2;
    }
    
    // check final parameters
    Eigen::Vector2d P1 = e1.point_at(theta1);
    Eigen::Vector2d P2 = e2.point_at(theta2);

    // Sanity check: are the points outside of the other ellipse?
    if (e2.is_inside(P1) || e1.is_inside(P2)) return true;

    // Check if the distance^2 between points is less than h^2
    return (P1 - P2).squaredNorm() <= h_sq;
}

std::pair<double, double> EllipseBundle::starting_parameters(const Ellipse &e1, const Ellipse &e2) const {
    double theta1 = 0,theta2 = 0;
    double min_start_dist_sq = 1e10;

    //center test
    const Eigen::Vector2d diff = e2.center - e1.center;

    //parameter for e1
    Eigen::Vector2d dir1 = e1.parametrize_matrix.inverse() * diff;
    double t1_center = std::atan2(dir1.y(), dir1.x());

    //parameter for e2
    Eigen::Vector2d dir2 = e2.parametrize_matrix.inverse() * (-diff);
    double t2_center = std::atan2(dir2.y(), dir2.x());

    Eigen::Vector2d p1 = e1.point_at(t1_center);
    Eigen::Vector2d p2 = e2.point_at(t2_center);
    double dist_sq = (p1-p2).squaredNorm();

    if (dist_sq <= min_start_dist_sq) {
        min_start_dist_sq = dist_sq;
        theta1 = t1_center;
        theta2 = t2_center;
    }

    const std::vector<double> e1_adj_angles = adjacent_angles(t1_center);
    const std::vector<double> e2_adj_angles = adjacent_angles(t2_center);

    //check if best starting point changes
    for (const double t1 : e1_adj_angles){
        p1 = e1.point_at(t1);   
        for (const double t2 : e2_adj_angles){
            p2 = e2.point_at(t2);
            dist_sq = (p1-p2).squaredNorm();
            
            if (dist_sq <= min_start_dist_sq) {
                min_start_dist_sq = dist_sq;
                theta1 = t1;
                theta2 = t2;
            }
        }
    }

    return {theta1, theta2};
}

std::vector<double> EllipseBundle::adjacent_angles(const double theta) const {
    const std::vector<double> angles = {0.0, pi/2, pi, 3 * pi / 2};
    const double two_over_pi = 2.0 / pi;

    int index = static_cast<int>(theta * two_over_pi); //floors the result when casting
    index = index % 4; //range of atan is [-pi,pi]

    return {angles[index], angles[(index + 1) % 4]};
}
