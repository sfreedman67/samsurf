from collections import namedtuple
from functools import lru_cache

from sage.all import *

from bowman.triangle import Triangle
from bowman.halfplane import HalfPlane


class Hinge(namedtuple("Hinge", ["vectors", "id_edge", "id_edge_opp"])):
    __slots__ = ()

    @property
    def coordinates(self):
        return tuple(coord for vector in self.vectors for coord in vector)

    def __hash__(self):
        return hash((self.coordinates,))

    @classmethod
    def from_id_edge(cls, trin, id_edge):
        label_tri, label_edge = id_edge

        id_edge_opp = trin.gluings[id_edge]

        label_tri_opp, label_edge_opp = id_edge_opp

        tri = trin.triangles[label_tri]
        tri_opp = trin.triangles[label_tri_opp]

        edge = tri[label_edge]
        edge_opp = tri_opp[label_edge_opp]

        if edge != -edge_opp:
            raise ValueError("Edges either nonparallel or improperly oriented")

        v0 = tri[(label_edge + 1) % 3]
        v1 = edge_opp
        v2 = -tri_opp[(label_edge_opp - 1) % 3]

        return Hinge((v0, v1, v2), id_edge, id_edge_opp)

    @property
    def is_convex(self):
        v0, v1, v2 = self.vectors
        boundary = [v0, v1 - v0, v2 - v1, -v2]
        crosses = [w0x * w1y - w1x * w0y
                   for (w0x, w0y), (w1x, w1y)
                   in zip(boundary, boundary[1:] + boundary[:1])]

        all_positive = all(bool(cross > 0) for cross in crosses)
        all_negative = all(bool(cross < 0) for cross in crosses)

        return all_positive or all_negative

    def flip(self):
        v0, v1, v2 = self.vectors
        return Hinge((v1 - v0, v2 - v0, -v0), self.id_edge, self.id_edge_opp)

    @property
    def incircle_det(self):
        """(p2 is inside/on/outside oriented circle 0-P0-P1) iff (det </==/> 0) """
        return sage.all.matrix([[x, y, x ** 2 + y ** 2] for x, y in self.vectors]).determinant()

    @property
    @lru_cache(None)
    def _coefficients(self):
        (x0, y0), (x1, y1), (x2, y2) = self.vectors

        m02 = x1 * y2 - x2 * y1
        m12 = x0 * y2 - x2 * y0
        m22 = x0 * y1 - x1 * y0

        a = y0 ** 2 * m02 - y1 ** 2 * m12 + y2 ** 2 * m22
        b = 2 * (x0 * y0 * m02 - x1 * y1 * m12 + x2 * y2 * m22)
        c = x0 ** 2 * m02 - x1 ** 2 * m12 + x2 ** 2 * m22

        return a, b, c

    @property
    def halfplane(self):
        try:
            return HalfPlane.from_ineq(*self._coefficients)
        except ValueError:
            return None

    @property
    def triangle(self):
        v0, v1, v2 = self.vectors

        sides_ordered = sorted([(self.id_edge[1], -v1),
                                ((self.id_edge[1] + 1) % 3, v0),
                                ((self.id_edge[1] + 2) % 3, v1 - v0)])
        return Triangle(*(vector for _, vector in sides_ordered))

    @property
    def triangle_opp(self):
        v0, v1, v2 = self.vectors

        sides_ordered = sorted([(self.id_edge_opp[1], v1),
                                ((self.id_edge_opp[1] + 1) % 3, v2 - v1),
                                ((self.id_edge_opp[1] + 2) % 3, -v2)])

        return Triangle(*(vector for _, vector in sides_ordered))

    @property
    def ids_boundary(self):
        """return the edge IDs of the boundary of the hinge
        starting in the NE and moving Clockwise"""

        label_tri, label_edge = self.id_edge
        label_tri_opp, label_edge_opp = self.id_edge_opp

        SE = (label_tri, (label_edge + 1) % 3)
        NE = (label_tri, (label_edge + 2) % 3)
        NW = (label_tri_opp, (label_edge_opp + 1) % 3)
        SW = (label_tri_opp, (label_edge_opp + 2) % 3)

        return NE, SE, SW, NW

    def plot(self):
        v0, v1, v2 = self.vectors
        vertices_t1 = [sage.all.zero_vector(2), v0, v1]
        vertices_t2 = [sage.all.zero_vector(2), v1, v2]
        return sage.all.polygon2d(vertices_t1, fill=False).plot() + sage.all.polygon2d(vertices_t2, fill=False).plot()
