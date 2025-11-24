#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "ellipse.h"
#include "square_solver.h"

namespace py = pybind11;

PYBIND11_MODULE(square_solver, m) {
    m.doc() = "ManiFEM Ellipse Solver";

    py::class_<Ellipse>(m, "Ellipse")
        .def(py::init<double, double, double, double, double>(),
             py::arg("x"), py::arg("y"), py::arg("A"), py::arg("B"), py::arg("C"));

    py::class_<EllipseBundle>(m, "EllipseBundle")
        .def(py::init<>())
        .def("add", [](EllipseBundle &self, const Ellipse &e) {
            self.add(e);
        }, py::arg("ellipse"))
        .def("area", &EllipseBundle::area);

    m.def("objective", [](double heat_sources, double base_temp, const EllipseBundle& ellipses, bool mesh) {
        maniFEM::Function h_func = heat_sources;
        maniFEM::Function b_func = base_temp;
        return objective(h_func, b_func, ellipses, mesh);
    }, 
    py::arg("heat_source"), 
    py::arg("base_temp"), 
    py::arg("ellipses"), 
    py::arg("export_mesh"));    
}