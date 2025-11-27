"""
Python bindings for SquareSolver and EllipseBundle
"""
from __future__ import annotations
import numpy
import typing
__all__ = ['DomainConfig', 'Ellipse', 'EllipseBundle', 'SquareSolver']
class DomainConfig:
    @typing.overload
    def __init__(self) -> None:
        """
        Default constructor
        """
    @typing.overload
    def __init__(self, x_m: float, y_m: float, MW: float, ME: float, h_size: float, n_ellipses: int) -> None:
        """
        Parameterized constructor
        """
    def n_segments(self, arg0: float) -> int:
        ...
    @property
    def ME_x(self) -> float:
        ...
    @property
    def MW_x(self) -> float:
        ...
    @property
    def h(self) -> float:
        ...
    @property
    def num_ellipses(self) -> int:
        ...
    @property
    def x_max(self) -> float:
        ...
    @property
    def y_max(self) -> float:
        ...
class Ellipse:
    center: numpy.ndarray[numpy.float64[2, 1]]
    def __init__(self, x_i: float, y_i: float, A_i: float, B_i: float, C_i: float) -> None:
        ...
    def area(self) -> float:
        ...
    def derivative_at(self, arg0: float) -> numpy.ndarray[numpy.float64[2, 1]]:
        ...
    def meshify(self, arg0: float) -> None:
        ...
    def point_at(self, arg0: float) -> numpy.ndarray[numpy.float64[2, 1]]:
        ...
    @property
    def height(self) -> float:
        ...
    @property
    def width(self) -> float:
        ...
class EllipseBundle:
    def __init__(self, arg0: DomainConfig) -> None:
        ...
    def add(self, arg0: Ellipse) -> None:
        """
        Add an ellipse to the bundle (throws if invalid)
        """
    def area(self) -> float:
        ...
    @property
    def cfg(self) -> DomainConfig:
        ...
class SquareSolver:
    def __init__(self, cfg: DomainConfig, heat_sources: float, base_temp: float, penalization: float, export_mesh: bool) -> None:
        ...
    def solve(self, bundle: EllipseBundle) -> float:
        """
        Solves the system and returns the objective value (GIL is released)
        """
