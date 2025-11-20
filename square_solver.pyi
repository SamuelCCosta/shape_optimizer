"""
ManiFEM Ellipse Solver
"""
from __future__ import annotations
__all__ = ['Ellipse', 'EllipseBundle', 'objective']
class Ellipse:
    def __init__(self, x: float, y: float, A: float, B: float, C: float) -> None:
        ...
class EllipseBundle:
    def __init__(self) -> None:
        ...
    def add(self, ellipse: Ellipse) -> None:
        ...
    def area(self) -> float:
        ...
def objective(heat_source: float, base_temp: float, ellipses: EllipseBundle, export_mesh: bool) -> float:
    ...
