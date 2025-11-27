#include "ellipse.h"

#include <cmath>

// A(x-xc)^2 + 2B(x-xc)(y-yc) + C(y-yc)^2 = 1
Ellipse::Ellipse(double x, double y, double A_, double B_, double C_) : A(A_), B(B_), C(C_),
    mesh(maniFEM::tag::non_existent) {
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

    Eigen::Matrix2d M_inv = M.inverse();
    Eigen::LLT<Eigen::Matrix2d> llt_inv(M_inv);
    transform = llt_inv.matrixL();
}

void Ellipse::meshify(const double h) const {
    WorkingManifold RR2WM;
    Manifold &RR2 = RR2WM.current;
    Function xy = RR2.coordinates();
    Function x = xy[0], y = xy[1];

    double xc = center[0], yc = center[1];
    Function implicit_eq = A * (x-xc) * (x-xc) + 2 * B * (x-xc) * (y-yc) + C * (y-yc) * (y-yc);
    Manifold representation = RR2.implicit(implicit_eq == 1);
    representation.set_as_working_manifold();

    Eigen::Vector2d starting_point = point_at(0);

    double x_0 = starting_point.x(), y_0 = starting_point.y();
    Eigen::Vector2d starting_derivative = derivative_at(0);

    //queremos o simétrico para a malha ter orientação contrária ao quadrado
    double x_dir = -starting_derivative.x(), y_dir = -starting_derivative.y();

    Cell start(tag::vertex, tag::of_coords, {x_0,y_0}); //extremidade equivalente a (a,0)
    std::vector<double> direction = {x_dir, y_dir}; //direção equivalente a (0,-1)
    mesh = Mesh::Build(tag::frontal).entire_manifold().start_at(start).towards(direction).desired_length(h);
}

bool EllipseBundle::intersects(const Ellipse &e1, const Ellipse &e2){
    double x1 = e1.center[0], y1 = e1.center[1];
    double x2 = e2.center[0], y2 = e2.center[1];

    double val2 = e2.A*(x1-x2)*(x1-x2) + 2*e2.B*(x1-x2)*(y1-y2) + e2.C*(y1-y2)*(y1-y2);
    if (val2 <= 1.0) return true; //e1 center inside e2
    double val1 = e1.A*(x1-x2)*(x1-x2) + 2*e1.B*(x1-x2)*(y1-y2) + e1.C*(y1-y2)*(y1-y2);
    if (val1 <= 1.0) return true; //e2 center inside e1

    //b = bottom, t = top, l = left, r = right
    double e1_l_x = x1 - e1.width; 
    double e1_b_y = y1 - e1.height;
    double e1_r_x = x1 + e1.width;
    double e1_t_y = y1 + e1.height;

    double e2_l_x = x2 - e2.width; 
    double e2_b_y = y2 - e2.height;
    double e2_r_x = x2 + e2.width;
    double e2_t_y = y2 + e2.height;

    if (e1_r_x + cfg.h < e2_l_x || //e1 à esquerda de e2
        e2_r_x + cfg.h < e1_l_x || //e1 à direita de e2
        e1_t_y + cfg.h < e2_b_y || //e1 abaixo de e2
        e2_t_y + cfg.h < e1_b_y)   //e1 acima de e2
    {return false;}
    
    return robust_intersect(e1,e2);
}

//a mudança no futuro
bool EllipseBundle::robust_intersect(const Ellipse& e1, const Ellipse& e2) const {
    const double h_sq = cfg.h*cfg.h;
    const double eps = 1e-09;
    const double max_iterations = 15;

    //transform columns are the semiaxis vectors
    auto u1 = e1.transform.col(0), v1 = e1.transform.col(1);
    auto u2 = e2.transform.col(0), v2 = e2.transform.col(1);

    double theta1 = 0.0;
    double theta2 = 0.0;
    double min_start_dist = 1e15;

    //Check tips
    const double angles[] = {0.0, pi/2, pi, 3 * pi/2};

    for(int i=0; i<4; ++i) {
        double t1 = angles[i];
        Eigen::Vector2d p1 = e1.center + u1 * std::cos(t1) + v1 * std::sin(t1);
        
        for(int j=0; j<4; ++j) {
            double t2 = angles[j];
            Eigen::Vector2d p2 = e2.center + u2 * std::cos(t2) + v2 * std::sin(t2);
            
            double d2 = (p1 - p2).squaredNorm();
            if(d2 < min_start_dist) {
                min_start_dist = d2;
                theta1 = t1;
                theta2 = t2;
            }
        }
    }

    // center test
    Eigen::Vector2d C_vec = e2.center - e1.center;
    double base1 = std::atan2(C_vec.dot(v1), C_vec.dot(u1));
    double base2 = std::atan2((-C_vec).dot(v2), (-C_vec).dot(u2));
    
    Eigen::Vector2d p1_c = e1.center + u1 * std::cos(base1) + v1 * std::sin(base1);
    Eigen::Vector2d p2_c = e2.center + u2 * std::cos(base2) + v2 * std::sin(base2);
    if ((p1_c - p2_c).squaredNorm() < min_start_dist) {
        theta1 = base1;
        theta2 = base2;
    }
    std::cout << "theta1: " << theta1 << " theta2: " << theta2 << std::endl;

    // P(t) = C + u cos(t) + v sin(t), P'(t) = -u sin(t) + v cos(t)
    // P''(t) = C - P(t) or P''(t) = -u cos(t) - v sin(t)
    for (int i = 0; i < max_iterations; ++i){
        double c1 = std::cos(theta1), s1 = std::sin(theta1);
        double c2 = std::cos(theta2), s2 = std::sin(theta2);

        //Check if dist^2 is smaller than h^2
        Eigen::Vector2d P1 = e1.center + u1 * c1 + v1 * s1;
        Eigen::Vector2d P2 = e2.center + u2 * c2 + v2 * s2;
        Eigen::Vector2d D = P2 - P1;
        double dist2 = D.squaredNorm();
        if (dist2 < h_sq) return true;
        
        Eigen::Vector2d Pp1 = -u1 * s1 + v1 * c1;
        Eigen::Vector2d Pp2 = -u2 * s2 + v2 * c2;

        // we want r1,r2 = 0
        double r1 = D.dot(Pp1), r2 = -D.dot(Pp2); 
        // If the vectors are perpendicular, we found our minima
        if (std::abs(r1) < eps && std::abs(r2) < eps) break;

        Eigen::Vector2d Ppp1 = -u1 * c1 - v1 * s1;
        Eigen::Vector2d Ppp2 = -u2 * c2 - v2 * s2;

        //calculate Hessian matrix of f_i(t1,t2)= D.dot(PPi)
        // H00 = P1' P1' + D P1''
        double H00 = Pp1.squaredNorm() - D.dot(Ppp1);
        // H01 = - P1' P2'
        double H01 = -Pp1.dot(Pp2);
        // H11 = - P2' P2' + D P2''
        double H11 = Pp2.squaredNorm() + D.dot(Ppp2);

        double det = H00 * H11 - H01 * H01;
        if (std::abs(det) < eps) return true; //Degenerate Hessian

        // Update parameters
        // Cramer's rule
        // [ H00 H01 ] [ dt1 ] = [ -r1 ]
        // [ H01 H11 ] [ dt2 ] = [ -r2 ]

        double dtheta1 = (H11 * r1 - H01 * r2) / det;
        double dtheta2 = (H00 * r2 - H01 * r1) / det;

        theta1 += dtheta1, theta2 += dtheta2;
    }

    double c1 = std::cos(theta1), s1 = std::sin(theta1);
    double c2 = std::cos(theta2), s2 = std::sin(theta2);

    Eigen::Vector2d P1 = e1.center + u1 * c1 + v1 * s1;
    Eigen::Vector2d P2 = e2.center + u2 * c2 + v2 * s2;

    std::cout << std::sqrt((P2 - P1).squaredNorm()) << std::endl;

    return (P2 - P1).squaredNorm() < h_sq;
}

std::pair<double, double> EllipseBundle::get_initial_thetas(const Ellipse& e1, const Ellipse& e2) const {
    Eigen::Vector2d center_diff = e2.center - e1.center;
    
    Eigen::Matrix2d A1_inv = e1.transform.inverse();
    Eigen::Matrix2d A2_inv = e2.transform.inverse();

    Eigen::Vector2d v_target1 = A1_inv * center_diff;
    double theta1 = std::atan2(v_target1.y(), v_target1.x());

    Eigen::Vector2d v_target2 = A2_inv * (-center_diff);
    double theta2 = std::atan2(v_target2.y(), v_target2.x());

    return {theta1, theta2};
}

const Mesh EllipseBundle::total_mesh() const {
    std::vector<Mesh> mesh_bundle;
    mesh_bundle.reserve(cfg.num_ellipses);

    for (auto& ellipse : bundle){
        mesh_bundle.push_back(ellipse.get_mesh(cfg.h));
    }

    return Mesh::Build(tag::join).meshes(mesh_bundle);
}