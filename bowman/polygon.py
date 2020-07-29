import sage.all
from sage.all import *

from collections import namedtuple
import itertools

from context import bowman
import bowman.radical
from bowman import radical
import bowman.halfplane
from bowman import halfplane


class Point(namedtuple("Point", ["u", "v2"])):
    def __new__(cls, u, v2):
        if not isinstance(u, radical.Radical):
            self = super(Point, cls).__new__(
                cls, radical.Radical(u, 0, 0), v2)
        else:
            self = super(Point, cls).__new__(cls, u, v2)

        return self

    @staticmethod
    def CCW(p1, p2, p3):
        if any((p1 == p2, p1 == p3, p2 == p3)):
            raise ValueError("Can only determine CCW for 3 distinct points")

        elif all((p1 != oo, p2 != oo, p3 != oo)):
            return (p1.u < p2.u < p3.u) or (p2.u < p3.u < p1.u) or (p3.u < p1.u < p2.u)

        elif p1 == oo:
            return p2.u < p3.u

        elif p2 == oo:
            return p3.u < p1.u

        return p1.u < p2.u


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

    def reverse(self):
        return Edge(self.halfplane.reorient(), self.end, self.start)

    def plot(self):
        if self.start == oo:
            coord_start = oo
        else:
            coord_start = CC(self.start.u.value, QQbar(self.start.v2).sqrt())

        if self.end == oo:
            coord_end = oo
        else:
            coord_end = CC(self.end.u.value, QQbar(self.end.v2).sqrt())

        return HyperbolicPlane().UHP().get_geodesic(coord_start, coord_end).plot(axes=True)


class Polygon(namedtuple("Polygon", ["edges"])):
    def __key(self):
        return tuple(sorted((edge.start for edge in self.edges)))

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
        
        # TODO: this line is always problematic
        [head_idx] = [idx for idx, edge in enumerate(chain)
                      if halfplane.contains_point_on_boundary(edge.start)]

        chain_rotated = [*chain[head_idx:], *chain[:head_idx]]

        edge_new = Edge(halfplane,
                        chain_rotated[-1].end,
                        chain_rotated[0].start)

        return Polygon([*chain_rotated, edge_new])

    def partition_edges(self, halfplane):
        test_starts = [halfplane.contains_point(start)
                       for _, start, _ in self.edges]

        test_endpoints = list(zip(test_starts,
                                  test_starts[1:] + test_starts[:1]))
        edges_new = []
        for edge, tested_endpoints in zip(self.edges, test_endpoints):
            components_new = halfplane.intersect_edge(edge, *tested_endpoints)
            edges_new.extend(component for component in components_new
                             if isinstance(component, Edge))
        return edges_new

    def intersect_with_halfplane(self, halfplane):
        edges_new = self.partition_edges(halfplane)

        if edges_new == self.edges:
            return self
        elif not edges_new:
            return None

        return Polygon._close_edge_chain(edges_new, halfplane)

    def plot(self):
        return sum(edge.plot() for edge in self.edges)
