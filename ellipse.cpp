#include "ellipse.h"

#include <cmath>

// A(x-xc)^2 + 2B(x-xc)(y-yc) + C(y-yc)^2 = 1
Ellipse::Ellipse(double x, double y, double A_, double B_, double C_) : A(A_), B(B_), C(C_) {
    this-> det = A*C - B*B;
    if (A < 0 || det < 0){ 
        throw std::invalid_argument("The matrix is not positive definite");
    }

    this-> M << A, B, B, C;
    this-> center << x, y;

    //bounding boxes half heights
    this-> height = std::sqrt(A/det);
    this-> width = std::sqrt(C/det);
    //https://www.geometrictools.com/Documentation/RobustIntersectionOfEllipses.pdf

    //bounds checking
    double horizontal_margin = h + width;
    double vertical_margin = h + height;
    std::cout << "x_min: " << horizontal_margin << " y_min: " << vertical_margin << std::endl;
    std::cout << "x_max: " << 1 - horizontal_margin << " y_max: " << 1 - vertical_margin << std::endl;
    std::cout << "x: " << x << " y: " << y << std::endl;
    if ((y < vertical_margin) || (x < horizontal_margin) ||
    (y > 1 - vertical_margin) || (x > 1 - horizontal_margin) ) {
        throw std::invalid_argument("Ellipse does not fit in the unit square.");
    }
}

void Ellipse::meshify() const {
    WorkingManifold RR2WM;
    Manifold &RR2 = RR2WM.current;
    Function xy = RR2.coordinates();
    Function x = xy[0], y = xy[1];

    double xc = center[0], yc = center[1];
    if (! representation.has_value()){
        Function implicit_eq = A * (x-xc) * (x-xc) + 2 * B * (x-xc) * (y-yc) + C * (y-yc) * (y-yc);
        representation.emplace(RR2.implicit(implicit_eq == 1));
    }
    representation->set_as_working_manifold();

    Eigen::Vector2d starting_point = point_at(0);

    double x_0 = starting_point.x(), y_0 = starting_point.y();
    Eigen::Vector2d starting_derivative = derivative_at(0);

    //queremos o simétrico para a malha ter orientação contrária ao quadrado
    double x_dir = -starting_derivative.x(), y_dir = -starting_derivative.y();

    Cell start(tag::vertex, tag::of_coords, {x_0,y_0}); //extremidade equivalente a (a,0)
    std::vector<double> direction = {x_dir, y_dir}; //direção equivalente a (0,-1)
    mesh.emplace(Mesh::Build(tag::frontal).entire_manifold().start_at(start).towards(direction).desired_length(h));
}

bool EllipseBundle::intersects(const Ellipse &e1, const Ellipse &e2){
    //b = bottom, t = top, l = left, r = right
    double x1 = e1.center[0], y1 = e1.center[1];
    double x2 = e2.center[0], y2 = e2.center[1];
    double e1_l_x = x1 - e1.width; 
    double e1_b_y = y1 - e1.height;
    double e1_r_x = x1 + e1.width;
    double e1_t_y = y1 + e1.height;

    double e2_l_x = x2 - e2.width; 
    double e2_b_y = y2 - e2.height;
    double e2_r_x = x2 + e2.width;
    double e2_t_y = y2 + e2.height;

    if (e1_r_x < e2_l_x || //e1 à esquerda de e2
        e2_r_x < e1_l_x || //e1 à direita de e2
        e1_t_y < e2_b_y || //e1 abaixo de e2
        e2_t_y < e1_b_y)   //e1 acima de e2
    {return false;}
    
    return robust_intersect(e1,e2);
}

//provavelmente será mudado no futuro por algo com mais performance (maybe not)
bool EllipseBundle::robust_intersect(const Ellipse& e1, const Ellipse& e2) const {
    auto [theta1, theta2] = get_initial_thetas(e1, e2); //ponto inicial
    const double learning_rate = 5;
    const int max_iterations = 200;
    const double tolerance_sq = 1e-10;

    for (int i = 0; i < max_iterations; ++i) {
        Eigen::Vector2d p1 = e1.point_at(theta1);
        Eigen::Vector2d p2 = e2.point_at(theta2);
        Eigen::Vector2d dist_vec = p1 - p2;
        double grad1 = 2 * dist_vec.dot(e1.derivative_at(theta1));
        double grad2 = -2 * dist_vec.dot(e2.derivative_at(theta2));
        if (grad1*grad1 + grad2*grad2 < tolerance_sq) {
            std::cout << "número de iterações: " << i << std::endl;
            break;
        }
        theta1 -= learning_rate * grad1;
        theta2 -= learning_rate * grad2;
    }
    double min_sq_dist = (e1.point_at(theta1) - e2.point_at(theta2)).squaredNorm();
    
    std::cout << "min dist^2: " << min_sq_dist << std::endl;
    std::cout << "h^2: " << h*h << std::endl;
    return min_sq_dist <= h * h + 1e-9;
}

std::pair<double, double> EllipseBundle::get_initial_thetas(const Ellipse& e1, const Ellipse& e2) const {
    Eigen::Vector2d center_diff = e2.center - e1.center;
    const auto& A1 = e1.get_transform_matrix();
    const auto& A2 = e2.get_transform_matrix();
    
    Eigen::Matrix2d A1_inv = A1.inverse();
    Eigen::Matrix2d A2_inv = A2.inverse();

    Eigen::Vector2d v_target1 = A1_inv * center_diff;
    double theta1 = std::atan2(v_target1.y(), v_target1.x());

    Eigen::Vector2d v_target2 = A2_inv * (-center_diff);
    double theta2 = std::atan2(v_target2.y(), v_target2.x());

    return {theta1, theta2};
}

Mesh EllipseBundle::total_mesh(){
    std::vector<Mesh> mesh_bundle;
    mesh_bundle.reserve(num_ellipses);

    for (auto& ellipse : bundle){
        mesh_bundle.push_back(ellipse.get_mesh());
    }

    return Mesh::Build(tag::join).meshes(mesh_bundle);
}