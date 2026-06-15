#include "ellipse.h"
#include <iostream>
#include <random>

Ellipse::Ellipse(double x, double y, double A, double B, double C) {
    this-> quadratic_form << A, B, B, C;
    this-> center << x, y;
    double det = quadratic_form.determinant();
    if (A <= 0 || det <= 0) {
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
    std::vector<double> direction = {start_derivative.x(), start_derivative.y()}; // Reverse orientation relative to the square boundary

    Mesh mesh = Mesh::Build(tag::frontal).entire_manifold().start_at(start).towards(direction).desired_length(h);

    RR2.set_as_working_manifold(); // TODO: implement RAII

    return mesh;
}

bool EllipseBundle::intersects(const Ellipse &e1, const Ellipse &e2) const {
    // Inexpensive: check if centers are inside eachother
    if (e2.evaluate_at(e1.center) <= 1.0) {return true;}
    if (e1.evaluate_at(e2.center) <= 1.0) {return true;}

    // Semi-Inexpensive: bounding boxes
    constexpr double three_halfs = 1.5;
    Eigen::Vector2d h_vec(three_halfs * h, three_halfs * h);

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
    constexpr bool ROBUST_INTERSECT_VERBOSITY = false;
    const double three_halfs_times_h_sq = 1.5 * h * h;

    if (ROBUST_INTERSECT_VERBOSITY) {
    std::cout << "target dist^2: " << three_halfs_times_h_sq << std::endl;
    }

    // stopping criterion
    const double eps = 1e-10;
    const int max_iterations = 50;
    
    //Damping coefficients
    double lambda = 0.0001;
    double mu = 10.0;

    // We are minimizing 1/2*||e1(t1) - e2(t2)||^2 where t_i in [0,2*pi]

    auto [theta1, theta2] = starting_parameters(e1, e2);
    if (ROBUST_INTERSECT_VERBOSITY) { 
        std::cout << "starting parameters: " << theta1 << " " << theta2 << std::endl;
    }

    for (int i = 0; i < max_iterations; i++) {
        Eigen::Vector2d P1 = e1.point_at(theta1);
        Eigen::Vector2d P2 = e2.point_at(theta2);

        // Early exit
        if (e2.is_inside(P1) || e1.is_inside(P2)) {return true;}

        // Check if the distance^2 between points is less than h^2
        Eigen::Vector2d difference = P1 - P2;
        if (ROBUST_INTERSECT_VERBOSITY) {
        std::cout << "dist^2: " << difference.squaredNorm() << std::endl;
        }
        // Early exit
        if (difference.squaredNorm() <= three_halfs_times_h_sq) {return true;}

        Eigen::Vector2d P_prime_1 = e1.derivative_at(theta1);
        Eigen::Vector2d P_prime_2 = e2.derivative_at(theta2);

        // jacobian transpose
        Eigen::Matrix2d jacobian;
        jacobian << P_prime_1, - P_prime_2;

        // Gauss-Newton approximate hessian
        Eigen::Matrix2d hessian_mat = jacobian.transpose() * jacobian; //TODO: make it so it doesn't square the condition number
        Eigen::Vector2d gradient = jacobian.transpose() * difference;

        // Stopping criteria: check if gradient is small
        if (gradient.squaredNorm() < eps) {break;}

        // Apply Levenberg-Marquardt damping and solve the linear system
        Eigen::Matrix2d hessian_LM = hessian_mat + lambda * Eigen::Matrix2d::Identity(); //could be hessian_mat + lambda * diag(hessian_mat)
        Eigen::Vector2d delta = hessian_LM.ldlt().solve(-gradient);

        // Stopping criteria: check if move is small
        if (delta.squaredNorm() < eps) {break;}

        double theta_test_1 = theta1 + delta(0);
        double theta_test_2 = theta2 + delta(1);

        Eigen::Vector2d difference_test = e1.point_at(theta_test_1) - e2.point_at(theta_test_2);
        if (difference_test.squaredNorm() < difference.squaredNorm()) { //tend towards Newton method
            theta1 = theta_test_1;
            theta2 = theta_test_2;
            lambda /= mu;
        } else { //tend towards gradient descent, reject last iteration
            lambda *= mu;
        }
    }

    // check final parameters
    Eigen::Vector2d P1 = e1.point_at(theta1);
    Eigen::Vector2d P2 = e2.point_at(theta2);

    // Sanity check
    if (e2.is_inside(P1) || e1.is_inside(P2)) {return true;}

    // Check if the distance^2 between points is less than h^2
    if (ROBUST_INTERSECT_VERBOSITY) {
        std::cout << "final dist^2: " << (P1-P2).squaredNorm() << std::endl;
    }
    return (P1 - P2).squaredNorm() <= three_halfs_times_h_sq;
}

std::pair<double, double> EllipseBundle::starting_parameters(const Ellipse &e1, const Ellipse &e2) const {
    //center test
    const Eigen::Vector2d diff = e2.center - e1.center;

    //parameter for e1
    Eigen::Vector2d dir1 = e1.parametrize_matrix.inverse() * diff;
    double t1_center = std::atan2(dir1.y(), dir1.x());
    double theta1 = t1_center;

    //parameter for e2
    Eigen::Vector2d dir2 = e2.parametrize_matrix.inverse() * (-diff);
    double t2_center = std::atan2(dir2.y(), dir2.x());
    double theta2 = t2_center;


    Eigen::Vector2d p1 = e1.point_at(t1_center);
    Eigen::Vector2d p2 = e2.point_at(t2_center);
    double dist_sq = (p1-p2).squaredNorm();
    double min_start_dist_sq = dist_sq;

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

    int index = static_cast<int>(std::floor(theta * two_over_pi));
    index = (index % 4 + 4) % 4; //range of atan is [-pi,pi]
    return {angles[index], angles[(index + 1) % 4], theta};
}

void EllipseBundle::generate_random(unsigned int seed, size_t max_attempts) {
    std::mt19937 gen;
    if (seed != 0) {
        gen.seed(seed);
    } else {
        std::random_device rd;
        gen.seed(rd());
    }

    std::uniform_real_distribution<> dist_x(h, x_max - h);
    std::uniform_real_distribution<> dist_y(h, y_max - h);
    // TODO: change how those work depend on domain
    std::uniform_real_distribution<> dist_AC(10.0, 300.0);
    std::uniform_real_distribution<> dist_B(-100.0, 100.0);

    size_t attempts = 0;
    while (bundle.size() < num_ellipses && attempts < max_attempts) {
        attempts++;
        double x = dist_x(gen);
        double y = dist_y(gen);
        double A = dist_AC(gen);
        double B = dist_B(gen);
        double C = dist_AC(gen);

        try {
            this->add(Ellipse(x, y, A, B, C));
        } catch (const std::invalid_argument&) {
            // Silently ignore it and let the loop try again
        }
    }
    
    if (bundle.size() < num_ellipses) {
        std::cerr << "Warning: Reached max attempts (" << max_attempts 
                  << ") before placing all ellipses. Total: " << bundle.size() << "\n";
    }
}

void EllipseBundle::fill() {
    // fills with vertical ellipses
    // all of them are 1.5 * h away from each other, as well as the outer boundary
    double small_axis = (x_max - 2 * h) / (num_ellipses * 2) - 0.75 * h;
    assert(small_axis > 0);
    double target_large_axis = y_max / 2 - 1.5 * h;

    // due to eccentricity cosntraints, we need to cap the larger axis
    double max_large_axis = 2.7206 * small_axis;
    double large_axis = std::min(target_large_axis, max_large_axis);
    double A = 1. / (small_axis * small_axis);
    double C = 1. / (large_axis * large_axis);

    // double eccentricity = std::sqrt(1.0 - (C / A));
    // std::cout << "eccentricity: " << eccentricity << std::endl;

    double x_center_start = x_max / (num_ellipses * 2) + h;
    double step = (x_max - 2 * h) / num_ellipses;
    for (size_t i = 0; i < num_ellipses; i++) {
        double x_center = x_center_start + i * step;
        this->add(Ellipse(x_center, 0.5, A, 0.0, C));
    }
}