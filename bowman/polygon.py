import sage.all
from sage.all import *

from collections import namedtuple
import itertools

from bowman.radical import Radical
import bowman.halfplane as halfplane


class Point(namedtuple('Point', ['u', 'v2'])):
    __slots__ = ()

    def __new__(cls, u, v2):
        self = super(Point, cls).__new__(cls, u, v2)

        if u != oo and not isinstance(u, Radical):
            return Point.__new__(cls, Radical(u, 0, 0), v2)

        return self

    def __repr__(self):
        return f"Point({self.u}, {self.v2})"

    @staticmethod
    def CCW(p1, p2, p3):
        if any((p1 == p2, p1 == p3, p2 == p3)):
            raise ValueError("Can only determine CCW for 3 distinct points")
        elif any(getattr(pt, "v2") != 0 for pt in (p1, p2, p3)):
            raise ValueError("Can only determine CCW for boundary points")

        if all(not pt.is_infinity for pt in (p1, p2, p3)):
            return (p1 < p2 < p3) or (p2 < p3 < p1) or (p3 < p1 < p2)

        elif p1.is_infinity:
            return p2 < p3

        elif p2.is_infinity:
            return p3 < p1

        return p1 < p2

    @property
    def is_infinity(self):
        return self.u == oo and self.v2 == 0

    @staticmethod
    def _plug_point_into_halfplane(point, plane):
        a, b, c = plane
        u, v2 = point
        A, B, C = u

        # a[(A + B sqrt(C))^2 + v2] + b (A + B sqrt(C)) + c >= 0
        # --> A1 + B1 sqrt(C) >= 0

        A1 = a * (A**2 + B**2 * C + v2) + b * A + c
        B1 = a * (2 * A * B) + b * B

        return Radical(A1, B1, C)


class Edge(namedtuple("Edge", ['halfplane', 'start', 'end'])):
    __slots__ = ()

    def __new__(cls, halfplane, start, end):
        if start == end:
            return start

        self = super(Edge, cls).__new__(cls, halfplane, start, end)

        return self

    def __repr__(self):
        Ideal_descriptor = "Ideal" if self.is_ideal else ""
        return Ideal_descriptor + f"Edge({self.start}->{self.end})"

    @property
    def is_ideal(self):
        return self.halfplane is None

    @property
    def endpoints(self):
        return (self.start, self.end)

    def plot(self):
        coords_start = CC(self.start.u.value, self.start.v2)
        coords_end = CC(self.end.u.value, self.end.v2)

        if self.is_ideal or isinstance(self.halfplane, halfplane.Line):
            return line((coords_start, coords_end)).plot()

        coords_start = CC(self.start.u.value, self.start.v2)
        coords_end = CC(self.end.u.value, self.end.v2)
        return HyperbolicPlane().UHP().get_geodesic(coords_start, coords_end).plot()


def intersect_polygon_halfplane(polygon, halfplane):
        intersection = itertools.chain.from_iterable(
            (halfplane.intersect_edge(edge) for edge in polygon))

        edges_new = [component for component in intersection if isinstance(component, Edge)]
        
        if edges_new == polygon:
            return polygon
        elif not edges_new:
            return None

        [head_idx]=[idx for idx, edge in enumerate(edges_new)
                      if halfplane.contains_point(edge.start, on_boundary=True)]

        edge_chain=[*edges_new[head_idx:], *edges_new[:head_idx]]

        edge_new=Edge(halfplane, edge_chain[-1].end, edge_chain[0].start)

        return [*edge_chain, edge_new]

def plot_polygon(polygon):
    return sum(edge.plot() for edge in polygon)
