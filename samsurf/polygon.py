import sage.all
from sage.all import *

from collections import namedtuple

from samsurf import radical
from samsurf import halfplane

class Point(namedtuple("Point", ["u", "v2"])):

    def __new__(cls, u, v2):
        if not isinstance(u, radical.Radical):
            self = super(Point, cls).__new__(
                cls, radical.Radical(u, QQ(0), QQ(0)), v2)
        else:
            self = super(Point, cls).__new__(cls, u, v2)

        return self

    def __repr__(self):
        return f"Point(u={QQbar(self.u.value)}, v2={QQbar(self.v2)})"

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

    def apply_mobius(self, m):
        from bowman import mobius
        return mobius.apply_mobius(m, self)


class Edge(namedtuple("Edge", ['halfplane', 'start', 'end'])):
    __slots__ = ()

    def __new__(cls, halfplane, start, end):
        if start == end:
            return start

        self = super(Edge, cls).__new__(cls, halfplane, start, end)

        return self

    def __repr__(self):
        ideal_descriptor = "Ideal" if self.is_ideal else ""
        return ideal_descriptor + f"Edge({self.start}->{self.end})"

    @classmethod
    def from_two_points(cls, p0, p1):
        return Edge(halfplane.HalfPlane.from_two_points(p0, p1), p0, p1)

    @property
    def is_ideal(self):
        return self.halfplane is None

    @property
    def endpoints(self):
        return self.start, self.end

    def reverse(self):
        return Edge(self.halfplane.reorient(), self.end, self.start)

    def angle(self, e2):
        """
        Compute the angle between two consecutive edges e1 and e2
        If their common endpoint is ideal, return 0
        :param self: Edge
        :param e2: Edge
        """
        if self.end != e2.start:
            raise ValueError("self and e2 are not consecutive")
        elif self.end == oo or self.end.v2 == 0:
            return 0
        px, py = RR(self.end.u.value), RR(self.end.v2).sqrt()

        v0x, v0y = e2.tangent_vector(px, py)
        v1x, v1y = self.reverse().tangent_vector(px, py)

        angle = RR(atan2(v1y * v0x - v0y * v1x, v0x * v1x + v0y * v1y))
        if angle < 0:
            angle += RR(2 * pi)
        return angle

    def tangent_vector(self, u, v):
        a, b, _ = self.halfplane
        if v == 0:
            raise ValueError("tangent vector not defined for endpoints")
        elif a == 0:
            return sage.all.vector([0, -b])
        else:
            a, b = RR(a), RR(b)
            num = -b - 2 * a * RR(u)
            denom = 2 * v
            dv_du = num / denom
            return sage.all.vector([a, dv_du])

    def apply_mobius(self, m):
        from bowman import mobius
        return Edge(self.halfplane.apply_mobius(m),
                    mobius.apply_mobius(m, self.start),
                    mobius.apply_mobius(m, self.end))

    def plot(self):
        return HyperbolicPlane().UHP().get_geodesic(*self.coordinates).plot(boundary=False, axes=True)

    @property
    def coordinates(self):
        coordinates = []
        for endpoint in [self.start, self.end]:
            if endpoint == oo:
                coordinate = oo
            else:
                x_val = 0 if QQbar(endpoint.u.value).is_zero() else QQbar(endpoint.u.value)
                y_val = 0 if QQbar(endpoint.v2).is_zero() else QQbar(sqrt(endpoint.v2))
                coordinate = CC(x_val, y_val)
            coordinates.append(coordinate)

        return coordinates


class ChainHasMultipleHeadsError(Exception):
    pass


class Polygon(namedtuple("Polygon", ["edges"])):

    def __key(self):
        return tuple(sorted(((edge.start == oo, edge.start) for edge in self.edges)))

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other):
        if isinstance(other, Polygon):
            return self.__key() == other.__key()
        return NotImplemented

    def __contains__(self, item):
        return item in self.edges

    def __iter__(self):
        return iter(self.edges)

    def __len__(self):
        return len(self.edges)

    @property
    def vertices(self):
        return tuple(edge.start for edge in self.edges)

    @staticmethod
    def _close_edge_chain(chain, halfplane):
        def is_closed(chain):
            return all(end == start
                       for (_, _, end), (_, start, _) in zip(chain, chain[1:] + chain[:1]))

        if is_closed(chain):
            return Polygon(chain)

        head_idxs = [idx for idx, edge in enumerate(chain)
                     if halfplane.contains_point_on_boundary(edge.start)]

        if len(head_idxs) != 1:
            heads = [chain[idx] for idx in head_idxs]
            raise ChainHasMultipleHeadsError(f"heads={heads}, chain={chain}")

        [head_idx] = head_idxs

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

    @property
    def area(self):
        angle_sum = sum(Edge.angle(e1, e2)
                        for e1, e2 in zip(self.edges, self.edges[1:] + self.edges[:1]))
        return RR(pi) * (len(self.edges) - 2) - angle_sum

    def plot(self):
        return sum(edge.plot() for edge in self.edges)

    def contains_point(self, point):
        return all(x.halfplane.contains_point(point) for x in self.edges)
