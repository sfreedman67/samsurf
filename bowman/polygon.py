import sage.all
from sage.all import *

from collections import namedtuple
import itertools

from context import bowman
import bowman.radical
from bowman import radical
import bowman.halfplane
from bowman import halfplane


class Point:
    def __init__(self, u, v2):
        if not isinstance(u, radical.Radical) and u != oo:
            self.u = radical.Radical(u, QQ(0), QQ(0))
            self.v2 = v2
        else:
            self.u = u
            self.v2 = v2

    def __repr__(self):
        return f"Point({self.u}, {self.v2})"

    def __iter__(self):
        return iter((self.u, self.v2))

    def __key(self):
        return (self.u, self.v2)

    def __hash__(self):
        return hash(self.__key())

    @staticmethod
    def CCW(p1, p2, p3):
        if any((p1 == p2, p1 == p3, p2 == p3)):
            raise ValueError("Can only determine CCW for 3 distinct points")
        elif any(getattr(pt, "v2") != QQ(0) for pt in (p1, p2, p3)):
            raise ValueError("Can only determine CCW for boundary points")

        if all(not pt.is_infinity for pt in (p1, p2, p3)):
            return (p1.u < p2.u < p3.u) or (p2.u < p3.u < p1.u) or (p3.u < p1.u < p2.u)

        elif p1.is_infinity:
            return p2.u < p3.u

        elif p2.is_infinity:
            return p3.u < p1.u

        return p1.u < p2.u

    @property
    def is_infinity(self):
        return not isinstance(self.u, radical.Radical) and self.u == oo and self.v2 == QQ(0)

    def __eq__(self, other):
        if self.is_infinity:
            return other.is_infinity
        elif other.is_infinity:
            return self.is_infinity
        return self.__key() == other.__key()


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

    def __hash__(self):
        return hash((self.start, self.end))

    @property
    def is_ideal(self):
        return self.halfplane is None

    @property
    def endpoints(self):
        return (self.start, self.end)

    def reverse(self):
        return Edge(self.halfplane.reverse(), self.end, self.start)

    @property
    def tangent_vector(self):
        plane = self.halfplane

        if plane is None:
            return vector([1, 0])

        elif isinstance(plane, halfplane.Circle):
            b = plane.b / plane.a

            u, v2 = self.start
            return vector([radical.Radical(0, -1, v2), u + b / 2])

        elif self.start.is_infinity:
            return vector([1, 0])
        else:
            return vector([0, 1])

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

    @property
    def vertices(self):
        return (edge.start for edge in self.edges)

    def __key(self):
        return tuple(sorted(self.vertices, key=lambda vertex: (vertex.is_infinity, vertex.u, vertex.v2)))

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other):
        if isinstance(other, Polygon):
            return self.__key() == other.__key()
        return NotImplemented

    @staticmethod
    def _close_edge_chain(chain, halfplane):
        def is_closed(chain):
            return all(snd == fst for (_, _, snd), (_, fst, _) in zip(chain, chain[1:] + chain[:1]))

        if is_closed(chain):
            return Polygon(chain)

        [head_idx] = [idx for idx, edge in enumerate(chain)
                      if halfplane.contains_point_on_boundary(edge.start)]

        chain_rotated = [*chain[head_idx:], *chain[:head_idx]]

        edge_new = Edge(halfplane,
                        chain_rotated[-1].end,
                        chain_rotated[0].start)

        return Polygon([*chain_rotated, edge_new])

    def intersect_with_halfplane(self, halfplane):
        edges_new = [component for edge in self.edges
                     for component in halfplane.intersect_edge(edge)
                     if isinstance(component, Edge)]

        if edges_new == self.edges:
            return self
        elif not edges_new:
            return None

        return Polygon._close_edge_chain(edges_new, halfplane)

    def plot(self):
        return sum(edge.plot() for edge in self.edges)
