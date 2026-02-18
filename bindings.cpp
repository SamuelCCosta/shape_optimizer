#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "ellipse.h"
#include "square_solver.h"

namespace py = pybind11;

PYBIND11_MODULE(square_solver, m) {
    m.doc() = "Bindings for Ellipse and SquareSolver using maniFEM";

    // --- Ellipse Bindings ---
    py::class_<Ellipse>(m, "Ellipse")
        .def(py::init<double, double, double, double, double>(),
             py::arg("x"), py::arg("y"), py::arg("A"), py::arg("B"), py::arg("C"))
        
        // Member variables
        .def_readwrite("center", &Ellipse::center)
        .def_readwrite("quadratic_form", &Ellipse::quadratic_form)
        .def_readwrite("parametrize_matrix", &Ellipse::parametrize_matrix)
        .def_readwrite("bounds", &Ellipse::bounds)

        // Methods
        .def("area", &Ellipse::area)
        .def("width", &Ellipse::width)
        .def("height", &Ellipse::height)
        .def("bounding_half", &Ellipse::bounding_half)
        .def("evaluate_at", &Ellipse::evaluate_at, py::arg("point"))
        .def("is_inside", &Ellipse::is_inside, py::arg("point"))

        // Specific Overloads (Exposing only the 'double' variants)
        .def("point_at", 
             py::overload_cast<const double>(&Ellipse::point_at, py::const_), 
             py::arg("t"))
        .def("derivative_at", 
             py::overload_cast<const double>(&Ellipse::derivative_at, py::const_), 
             py::arg("t"))
        .def("second_derivative_at", 
             py::overload_cast<const double>(&Ellipse::second_derivative_at, py::const_), 
             py::arg("t"));


    // --- EllipseBundle Bindings ---
    py::class_<EllipseBundle>(m, "EllipseBundle")
        .def(py::init<std::map<std::string, double>&, const double, const size_t>(),
             py::arg("geometric_config"), py::arg("h"), py::arg("num_ellipses"))
        
        // Member variables
        .def_readwrite("bundle", &EllipseBundle::bundle)
        .def_readonly("x_max", &EllipseBundle::x_max)
        .def_readonly("y_max", &EllipseBundle::y_max)
        .def_readonly("MW_x", &EllipseBundle::MW_x)
        .def_readonly("ME_x", &EllipseBundle::ME_x)
        .def_readonly("h", &EllipseBundle::h)

        // Methods
        .def("add", &EllipseBundle::add, py::arg("new_ellipse"))
        .def("area", &EllipseBundle::area);


    // --- SquareSolver Bindings ---
    py::class_<SquareSolver>(m, "SquareSolver")
        // Binding the constructor that uses doubles for boundary conditions
        .def(py::init<std::map<std::string, double>&, const double, const double, const double, const double, const bool, const bool>(),
             py::arg("geometric_config"), 
             py::arg("h"),
             py::arg("heat_sources"), 
             py::arg("base_temp"),
             py::arg("penalization"), 
             py::arg("export_domain"), 
             py::arg("export_mesh"))
        
        // Member variables (subset based on public visibility)
        .def_readwrite("h", &SquareSolver::h)
        .def_readwrite("heat_sources", &SquareSolver::heat_sources)
        .def_readwrite("base_temp", &SquareSolver::base_temp)
        .def_readwrite("penalization", &SquareSolver::penalization)
        .def_readonly("export_domain", &SquareSolver::export_domain)
        .def_readonly("export_result", &SquareSolver::export_result)

        // Methods
        .def("solve", &SquareSolver::solve, py::arg("bundle"));
}