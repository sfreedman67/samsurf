import sage.all
from sage.all import *

from collections import namedtuple
import itertools

from context import bowman 

from bowman.radical import Radical
import bowman.halfplane as halfplane


class Point:
    def __init__(self, u, v2):
        if not isinstance(u, Radical) and u != oo:
            self.u = Radical(u, QQ(0), QQ(0))
            self.v2 = v2
        else:
            self.u = u
            self.v2 = v2

    def __repr__(self):
        return f"Point({self.u}, {self.v2})"

    def __iter__(self):
        return iter((self.u, self.v2))

    @staticmethod
    def CCW(p1, p2, p3):
        if any((p1 == p2, p1 == p3, p2 == p3)):
            raise ValueError("Can only determine CCW for 3 distinct points")
        elif any(getattr(pt, "v2") != QQ(0) for pt in (p1, p2, p3)):
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
        return not isinstance(self.u, Radical) and self.u == oo and self.v2 == QQ(0)

    def __eq__(self, other):
        if self.is_infinity:
            return other.is_infinity
        elif other.is_infinity:
            return self.is_infinity
        return (self.u, self.v2) == (other.u, other.v2)

    def __lt__(self, other):
        if self.is_infinity or other.is_infinity:
            raise ValueError("Can't compare infinity")
        return (self.u, self.v2) < (other.u, other.v2)


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
        if self.start.is_infinity:
            coord_start = oo
        else:
            coord_start = CC(self.start.u.value, QQbar(self.start.v2).sqrt())

        if self.end.is_infinity:
            coord_end = oo
        else:
            coord_end = CC(self.end.u.value, QQbar(self.end.v2).sqrt())

        return HyperbolicPlane().UHP().get_geodesic(coord_start, coord_end).plot(axes=True)

class Polygon(namedtuple("Polygon", ["edges"])):

    # TODO: is_nontriv property

    # TODO: two polys eq if same edges
    
    def intersect_with_halfplane(self, halfplane):
        intersection = itertools.chain.from_iterable(
            (halfplane.intersect_edge(edge) for edge in self.edges))

        edges_new = [component for component in intersection
                     if isinstance(component, Edge)]

        if edges_new == self.edges:
            return self
        elif not edges_new:
            return None

        [head_idx] = [idx for idx, edge in enumerate(edges_new)
                      if halfplane.contains_point_on_boundary(edge.start)]

        edge_chain = [*edges_new[head_idx:], *edges_new[:head_idx]]

        if edge_chain[0].start == edge_chain[-1].end:
            return Polygon(edge_chain)

        edge_new = Edge(halfplane, edge_chain[-1].end, edge_chain[0].start)

        return Polygon([*edge_chain, edge_new])


    def plot_polygon(self):
        return sum(edge.plot() for edge in self.edges)
