#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "domain_config.h"
#include "ellipse.h"
#include "square_solver.h"

namespace py = pybind11;

PYBIND11_MODULE(square_solver, m) {
    m.doc() = "Python bindings for SquareSolver and EllipseBundle";

    // --- DomainConfig Binding ---
    py::class_<DomainConfig>(m, "DomainConfig")
        .def(py::init<>(), "Default constructor")
        .def(py::init<double, double, double, double, double, size_t>(),
             py::arg("x_m"), py::arg("y_m"), py::arg("MW"), py::arg("ME"), 
             py::arg("h_size"), py::arg("n_ellipses"),
             "Parameterized constructor")
        .def_readonly("x_max", &DomainConfig::x_max)
        .def_readonly("y_max", &DomainConfig::y_max)
        .def_readonly("MW_x", &DomainConfig::MW_x)
        .def_readonly("ME_x", &DomainConfig::ME_x)
        .def_readonly("h", &DomainConfig::h)
        .def_readonly("num_ellipses", &DomainConfig::num_ellipses)
        .def("n_segments", &DomainConfig::n_segments);

    // --- Ellipse Binding ---
    py::class_<Ellipse>(m, "Ellipse")
        .def(py::init<double, double, double, double, double>(),
             py::arg("x_i"), py::arg("y_i"), py::arg("A_i"), py::arg("B_i"), py::arg("C_i"))
        .def("area", &Ellipse::area)
        .def("point_at", &Ellipse::point_at)
        .def("derivative_at", &Ellipse::derivative_at)
        .def_readwrite("center", &Ellipse::center) 
        .def_readonly("height", &Ellipse::height)
        .def_readonly("width", &Ellipse::width);

    // --- EllipseBundle Binding ---
    py::class_<EllipseBundle>(m, "EllipseBundle")
        .def(py::init<const DomainConfig&>())
        .def("add", [](EllipseBundle &self, const Ellipse &e) {
            self.add(e);
        }, "Add an ellipse to the bundle (throws if invalid)")
        .def("area", &EllipseBundle::area)
        .def_readonly("cfg", &EllipseBundle::cfg);

    // --- SquareSolver Binding ---
    py::class_<SquareSolver>(m, "SquareSolver")
        .def(py::init<DomainConfig, double, double, double, bool>(),
             py::arg("cfg"), 
             py::arg("heat_sources"), 
             py::arg("base_temp"), 
             py::arg("penalization"), 
             py::arg("export_mesh"))
        .def("solve", &SquareSolver::solve, 
             py::arg("bundle"), 
             py::call_guard<py::gil_scoped_release>(),
             "Solves the system and returns the objective value (GIL is released)");
}