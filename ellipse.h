#ifndef ELLIPSE_H
#define ELLIPSE_H

#include "maniFEM.h"
#include "constants.h"
#include <optional>
#include <numbers>

class WorkingManifold{
    public:
        WorkingManifold() = default;
        ~WorkingManifold() {current.set_as_working_manifold();}

        WorkingManifold(const WorkingManifold&) = delete;
        WorkingManifold& operator=(const WorkingManifold&) = delete;

        Manifold &current = Manifold::working;
};

class Ellipse {
    public:
        double xc, yc, a, b, theta; 
        std::optional<Manifold> representation;
        std::optional<Mesh> mesh;

        Ellipse(double xc_i, double yc_i, double a_i, double b_i, double theta_i);

        void meshify();

        Mesh get_mesh() {
            if (! mesh.has_value()){
                meshify();
            }
            return mesh.value();
        }
        
    private:
        void check_bounds(double xc_i, double yc_i, double a_i, double b_i, double theta_i) const;
};
#endif