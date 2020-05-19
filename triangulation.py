import sage.all
from sage.all import *

import flatsurf as fs

from itertools import product

import unittest


class Triangle:

    def __init__(self, edges):
        if sum(edges) != vector([0, 0]):
            raise ValueError("sides are oriented incorrectly")
        else:
            self.edges = edges

    def __repr__(self):
        return f'Triangle{self.edges[0], self.edges[1], self.edges[2]}'

    def __eq__(self, other):
        if isinstance(other, Triangle):
            return self.edges == other.edges
        return NotImplemented

    def apply_matrix(self, M):
        transformed_edges = [M * edge for edge in self.edges]
        return Triangle(transformed_edges)

class Hinge:

    def __init__(self, v1, v2, v3):
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3


class Triangulation:

    def __init__(self, triangles, gluings, base_ring=None, generator=None):
        self.triangles = triangles
        self.gluings = gluings
        self.ring = base_ring
        self.gen = generator

    def __repr__(self):
        return f'Triangulation(edges={self.triangles}, gluings={self.gluings})'

    def __eq__(self, other):
        if isinstance(other, Triangulation):
            return self.triangles == other.triangles and self.gluings == other.gluings
        return NotImplemented

    @classmethod
    def from_flatsurf(cls, X):
        DT = X.delaunay_triangulation()

        DT_polygons = [DT.polygon(i) for i in range(DT.num_polygons())]

        triangles = [Triangle([vector(edge) for edge in polygon.edges()])
                     for polygon in DT_polygons]

        gluings = {edge[0]: edge[1] for edge in DT.edge_iterator(gluings=True)}

        ring = DT.base_ring()

        gen = ring.gen()

        return Triangulation(triangles, gluings, ring, gen)

    @classmethod
    def square_torus(cls):
        return cls.from_flatsurf(fs.translation_surfaces.square_torus())

    @classmethod
    def regular_octagon(cls):
        return cls.from_flatsurf(fs.translation_surfaces.regular_octagon())

    @classmethod
    def arnoux_yoccoz(cls, g):
        if g < 3:
            raise ValueError("g must be >= 3")
        return cls.from_flatsurf(fs.translation_surfaces.arnoux_yoccoz(g))

    @classmethod
    def octagon_and_squares(cls):
        return cls.from_flatsurf(fs.translation_surfaces.octagon_and_squares())

    def num_triangles(self):
        return len(self.triangles)

    def opposite_edge(self, tri_lab, edge_lab):
        return self.gluings[(tri_lab, edge_lab)]

    def edges(self, gluings=False):
        edges = []

        for edge in product(range(len(self.triangles)), range(3)):
            opposite_edge = self.opposite_edge(edge[0], edge[1])
            if not (edge in edges or opposite_edge in edges):
                edges.append(edge)

        if gluings:
            edges = [(edge, self.opposite_edge(edge[0], edge[1]))
                     for edge in edges]
        return edges

    def apply_matrix(self, M):
        sheared_triangles = [triangle.apply_matrix(
            M) for triangle in self.triangles]
        return Triangulation(sheared_triangles, self.gluings)

    # TODO: Should hinge() remember its edges?
    def hinge(self, label_t1, label_e1):
        label_t2, label_e2 = self.opposite_edge(label_t1, label_e1)

        # by convention, the "first" triangle is the one with lower index
        if label_t1 > label_t2:
            label_t1, label_t2 = label_t2, label_t1
            label_e1, label_e2 = label_e2, label_e1

        triangle1, triangle2 = self.triangles[
            label_t1], self.triangles[label_t2]
        # edges are vectors
        edge1, edge2 = triangle1.edges[label_e1], triangle2.edges[label_e2]

        if edge1 != -edge2:
            raise ValueError(
                "Edges are either nonparallel or oriented incorrectly")

        v1 = triangle2.edges[(label_e2 + 1) % 3]
        v2 = edge1
        v3 = -triangle1.edges[(label_e1 - 1) % 3]

        if matrix([v2 - v1, v3 - v2]).determinant() < 0:
            v1, v3 = v3, v1

        return (v1, v2, v3)

    def hinges(self):
        return [self.hinge(label_t, label_e) for label_t, label_e in self.edges()]

    @staticmethod
    def edge_inequality(hinge):
        if matrix([hinge[1] - hinge[0], hinge[2] - hinge[1]]).determinant() < 0:
            raise ValueError("hinge is not oriented correctly")

        # TODO mess of a comment!
        # want det of matrix w/ row [xi + uyi, vyi, (xi + uyi)^2 + (vyi)^2 ]
        # expand out 3rd column, remove v from second column, use
        # multilinearity

        a = matrix([[v[0], v[1], v[1]**2] for v in hinge]).determinant()
        b = 2 * matrix([[v[0], v[1], v[0] * v[1]]
                        for v in hinge]).determinant()
        c = matrix([[v[0], v[1], v[0]**2] for v in hinge]).determinant()

        # TODO check degeneracy elsewhere...
        if a != 0:
            if b**2 - 4 * a * c <= 0:
                return None

        elif b == 0:
            if c >= 0:
                return None
            else:
                raise ValueError("a==b==0 and c < 0 is a degenerate inequality")

        return (a, b, c)

    def edge_inequalities(self):
        return [Triangulation.edge_inequality(hinge) for hinge in self.hinges() if Triangulation.edge_inequality(hinge) is not None]

    def is_non_degenerate(self):
        """
        TESTS::

            sage: Triangulation.is_non_degenerate(Triangulation.regular_octagon())
            False

        """
        return all(matrix([[v[0], v[1], v[0]**2 + v[1]**2] for v in hinge]).determinant() != 0 for hinge in self.hinges())

    @staticmethod
    def inequality_to_geodesic(ineq):
        a, b, c = ineq
        start, end = 0, 0
        if a == 0:
            if b == 0:
                raise ValueError("invalid inequality with a == b == 0")
            else:
                if b < 0:
                    start, end = -c / b, oo
                else:
                    start, end = oo, -c / b
        else:
            d = b**2 - 4 * a * c
            if d <= 0:
                raise ValueError("discriminant was non-positive")
            else:
                center = -b / (2 * a)
                left_root = center - (sqrt(d) / (2 * a))
                right_root = center + (sqrt(d) / (2 * a))
                if a * center**2 + b * center + c > 0:
                    start, end = right_root, left_root
                else:
                    start, end = left_root, right_root

        return HyperbolicPlane().UHP().get_geodesic(start, end)

    @staticmethod
    def edge_geodesic(hinge):
        return inequality_to_geodesic(edge_inequality(hinge))

    def edge_geodesics(self):
        return [inequality_to_geodesic(ineq) for ineq in self.edge_inequalities()]
