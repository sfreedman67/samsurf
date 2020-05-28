import sage.all
from sage.all import *

import flatsurf as fs

import itertools
import operator

import halfplane

import unittest


class Triangle:

    def __init__(self, edges):
        if sum(edges) != vector([0, 0]):
            raise ValueError("sides do not close up")
        elif bool(matrix([edges[0], edges[0] + edges[1]]).determinant() <= 0):
            raise ValueError("sides are not oriented correctly")
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

    def __init__(self, v0, v1, v2):
        self.vectors = [v0, v1, v2]

    def __repr__(self):
        return f"Hinge{self.vectors[0], self.vectors[1], self.vectors[2]}"

    def __eq__(self, other):
        if isinstance(other, Hinge):
            return self.vectors == other.vectors
        return NotImplemented

    @classmethod
    def from_triangles(cls, triangle1, label_e1, triangle2, label_e2):
        edge1, edge2 = triangle1.edges[label_e1], triangle2.edges[label_e2]

        if edge1 != -edge2:
            raise ValueError(
                "Edges are either nonparallel or oriented incorrectly")

        v1 = triangle2.edges[(label_e2 + 1) % 3]
        v2 = edge1
        v3 = -triangle1.edges[(label_e1 - 1) % 3]

        return Hinge(v1, v2, v3)

    def incircle_test(self):
        return matrix([[v[0], v[1], v[0]**2 + v[1]**2] for v in self.vectors]).determinant()

    def is_convex(self):
        pass

    def edge_inequality(self):
        a = matrix([[v[0], v[1], v[1]**2]
                    for v in self.vectors]).determinant()
        b = 2 * matrix([[v[0], v[1], v[0] * v[1]]
                        for v in self.vectors]).determinant()
        c = matrix([[v[0], v[1], v[0]**2]
                    for v in self.vectors]).determinant()

        return (a, b, c)

    def edge_geodesic(self):
        return halfplane.inequality_to_geodesic(self.edge_inequality())


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

        for edge in itertools.product(range(len(self.triangles)), range(3)):
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
        return Triangulation(sheared_triangles, self.gluings, self.ring, self.gen)

    def hinge(self, label_t1, label_e1):
        label_t2, label_e2 = self.opposite_edge(label_t1, label_e1)
        return Hinge.from_triangles(self.triangles[label_t1], label_e1, self.triangles[label_t2], label_e2)

    def hinges(self):
        return [self.hinge(*edge) for edge in self.edges()]

    # TODO: make edge inequalites and edge_geodesics into one object, a halfplane

    def edge_inequalities(self):
        return list(set(hinge.edge_inequality() for hinge in self.hinges()
                        if not halfplane.is_trivial(hinge.edge_inequality())))

    def is_delaunay(self, strict=False):
        compare = operator.ge if not strict else operator.gt
        return all(bool(compare(hinge.incircle_test(), 0)) for hinge in self.hinges())

    def edge_geodesics(self):
        return [halfplane.inequality_to_geodesic(*ineq) for ineq in self.edge_inequalities()]

    def plot_geodesics(self, count=None):
        return sum(itertools.islice(map(halfplane.plot, self.edge_geodesics()), count))
