from sage.all import *
import flatsurf as fs
from itertools import product
import unittest


class Triangle:

    def __init__(self, edges):
        if edges[0] + edges[1] != -edges[2]:
            raise ValueError("sides are oriented incorrectly")
        else:
            self.edges = edges

    def edge(self, n):
        return self.edges[n]

    def __repr__(self):
        return f'Triangle{self.edge(0), self.edge(1), self.edge(2)}'


class Triangulation:

    def __init__(self, triangles, gluings):
        self.triangles = triangles
        self.gluings = gluings

    def __repr__(self):
        return f'Triangulation(edges={self.triangles}, gluings={self.gluings})'

    def triangle(self, n):
        return self.triangles[n]

    def opposite_edge(self, edge):
        return self.gluings[edge]

    def edge_reps(self):
        edge_reps = set()

        for edge in product(range(len(self.triangles)), range(3)):
            opposite_edge = self.opposite_edge(edge)
            if not (edge in edge_reps or opposite_edge in edge_reps):
                edge_reps.add(edge)

        return edge_reps

    def hinge(self, edge):
        e1, e2 = edge, self.opposite_edge(edge)

        t_lab1, t_lab2 = e1[0], e2[0]
        e_lab1, e_lab2 = e1[1], e2[1]
        tri1, tri2 = self.triangle(t_lab1), self.triangle(t_lab2)
        # edges are vectors
        e1, e2 = tri1.edge(e_lab1), tri2.edge(e_lab2)

        if e1 != -e2:
            raise ValueError(
                "Edges are either nonparallel or oriented incorrectly")

        v2 = e2
        v1 = tri2.edge((e_lab2 + 1) % 3)
        v3 = tri1.edge((e_lab1 - 1) % 3)

        # want p1 -> p2, p2 -> p3 to be an oriented basis
        if matrix([v_2 - v_1, v3 - v2]).determinant() < 0:
            v1, v3 = v3, v1

        return (v1, v2, v3)

    def hinges(self):
        return {self.hinge(self, edge) for edge in self.edge_reps()}

    @classmethod
    def from_flatsurf(cls, X):
        DT = X.delaunay_triangulation()

        DT_polygons = [DT.polygon(i) for i in range(DT.num_polygons())]

        triangles = [Triangle([vector(edge) for edge in polygon.edges()])
                     for polygon in DT_polygons]

        gluings = {edge[0]: edge[1] for edge in DT.edge_iterator(gluings=True)}

        return Triangulation(triangles, gluings)


class TriangulationTests(unittest.TestCase):

    def testEdgeReps(self):
        T1 = Triangulation.from_flatsurf(fs.translation_surfaces.square_torus())
        self.assertEqual(T1.edge_reps(), {(0, 0), (0, 1), (0, 2)})

        T2 = Triangulation.from_flatsurf(
            fs.translation_surfaces.regular_octagon())
        print(T2)
        self.assertEqual(T2.edge_reps(), {(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (3, 0), (4, 0)})


if __name__ == "__main__":
    unittest.main(verbosity=2)
