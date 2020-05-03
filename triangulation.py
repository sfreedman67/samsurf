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
        return f'Triangle{self.edge(0), self.edge(1), self.edge(2)}'

    def __eq__(self, other):
        if isinstance(other, Triangle):
            return self.edges == other.edges
        return NotImplemented

    def edge(self, n):
        return self.edges[n]

    def apply_matrix(self, M):
        transformed_edges = [M * edge for edge in self.edges]
        return Triangle(transformed_edges)


class Triangulation:

    def __init__(self, triangles, gluings, generator=None):
        self.triangles = triangles
        self.gluings = gluings
        if gen is not None:
            self.gen = generator

    def __repr__(self):
        return f'Triangulation(edges={self.triangles}, gluings={self.gluings})'

    def __eq__(self, other):
        if isinstance(other, Triangulation):
            return self.triangles == other.triangles and self.gluings == other.gluings
        return NotImplemented

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

    def hinge(self, edge_data):
        opp_edge_data = self.opposite_edge(edge_data)

        t_lab1, t_lab2 = edge_data[0], opp_edge_data[0]
        e_lab1, e_lab2 = edge_data[1], opp_edge_data[1]
        tri1, tri2 = self.triangle(t_lab1), self.triangle(t_lab2)
        # edges are vectors
        e1, e2 = tri1.edge(e_lab1), tri2.edge(e_lab2)

        if e1 != -e2:
            raise ValueError(
                "Edges are either nonparallel or oriented incorrectly")

        # the hinge is "centered" at e1, meaning that e1 is the central vector
        v1 = tri2.edge((e_lab2 + 1) % 3)
        v2 = e1
        v3 = -tri1.edge((e_lab1 - 1) % 3)

        # want p1 -> p2, p2 -> p3 to be an oriented basis
        if matrix([v2 - v1, v3 - v2]).determinant() < 0:
            v1, v3 = v3, v1

        return (v1, v2, v3)

    @staticmethod
    def edge_inequality(hinge):
        if matrix([hinge[1] - hinge[0], hinge[2] - hinge[1]]).determinant() < 0:
            raise ValueError("hinge is not oriented correctly")

        # want det of matrix w/ row [xi + uyi, vyi, (xi + uyi)^2 + (vyi)^2 ]
        # expand out 3rd column, remove v from second column, use
        # multilinearity

        a = matrix([[v[0], v[1], v[1]**2] for v in hinge]).determinant()
        b = 2 * matrix([[v[0], v[1], v[0] * v[1]]
                        for v in hinge]).determinant()
        c = matrix([[v[0], v[1], v[0]**2] for v in hinge]).determinant()

        # return an inequality a(u^2 + v^2) + bu + c >= 0
        return (a, b, c)

    def edge_inequalities(self):
        return [Triangulation.edge_inequality(hinge) for hinge in self.hinges()]

    def hinges(self):
        return [self.hinge(edge) for edge in self.edge_reps()]

    def is_non_degenerate(self):
        return all(matrix([[v[0], v[1], v[0]**2 + v[1]**2] for v in hinge]).determinant() != 0 for hinge in self.hinges())

    def apply_matrix(self, M):
        sheared_triangles = [triangle.apply_matrix(
            M) for triangle in self.triangles]
        return Triangulation(sheared_triangles, self.gluings)

    @classmethod
    def from_flatsurf(cls, X):
        DT = X.delaunay_triangulation()

        DT_polygons = [DT.polygon(i) for i in range(DT.num_polygons())]

        triangles = [Triangle([vector(edge) for edge in polygon.edges()])
                     for polygon in DT_polygons]

        gluings = {edge[0]: edge[1] for edge in DT.edge_iterator(gluings=True)}

        gen = DT.base_ring().gen()

        return Triangulation(triangles, gluings, gen)

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


class TestApplyMatrixToTriangle(unittest.TestCase):

    def setUp(self):
        self.torus = Triangulation.square_torus()
        self.octagon = Triangulation.regular_octagon()
        self.a = self.octagon.gen
        self.A2 = matrix([[1, 2 * (self.a + 1)], [0, 1]])

    def test_shear_triangle0_in_torus(self):
        sq_t = Triangulation.square_torus()
        tri0 = sq_t.triangle(0)
        M = matrix([[2, 1], [1, 1]])
        sheared_triangle = Triangle(
            [vector([3, 2]), vector([-2, -1]), vector([-1, -1])])
        self.assertEqual(tri0.apply_matrix(M), sheared_triangle)

    def test_shear_triangle5_in_octagon(self):
        tri5 = self.octagon.triangle(5)
        w1 = vector([1 / 2 * (-8 - 5 * self.a), 1 / 2 * (-2 - self.a)])
        w2 = vector([6 + 4 * self.a, 1 + self.a])
        w3 = vector([1 / 2 * (-4 - 3 * self.a), -1 / self.a])
        sheared_tri5 = Triangle([w1, w2, w3])
        self.assertEqual(tri5.apply_matrix(self.A2), sheared_tri5)


class TestEdgeReps(unittest.TestCase):

    def test_regular_polygons(self):
        T1 = Triangulation.square_torus()

        self.assertEqual(T1.edge_reps(), {(0, 0), (0, 1), (0, 2)})

        T2 = Triangulation.regular_octagon()
        self.assertEqual(T2.edge_reps(), {
                         (0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (3, 0), (4, 0)})


class TestHinge(unittest.TestCase):

    def test_octagon_squares(self):
        T1 = Triangulation.octagon_and_squares()
        self.assertEqual(T1.hinge((1, 0)), (vector(
            [-2, 0]), vector([-2, -2]), vector([0, -2])))
        # TODO: write another test case


class TestEdgeInequality(unittest.TestCase):

    def test_incorrectly_oriented_hinge(self):
        h = (vector([0, 1]), vector([1, 1]), vector([1, 0]))
        self.assertRaises(ValueError, Triangulation.edge_inequality, h)

    def test_normal_hinge(self):
        h = (vector([2, 2]), vector([2, 4]), vector([1, 4]))
        self.assertEqual(Triangulation.edge_inequality(h), (-16, -32, -4))


class TestApplyMatrixToTriangulation(unittest.TestCase):

    def test_shear_torus(self):
        sq_t = Triangulation.square_torus()

        sheared_triangle0 = Triangle(
            [vector([3, 2]), vector([-2, -1]), vector([-1, -1])])
        sheared_triangle1 = Triangle(
            [vector([-3, -2]), vector([2, 1]), vector([1, 1])])
        sheared_triangulation = Triangulation(
            [sheared_triangle0, sheared_triangle1], sq_t.gluings)

        M = matrix([[2, 1], [1, 1]])
        self.assertEqual(sq_t.apply_matrix(M), sheared_triangulation)


class TestIsNonDegenerate(unittest.TestCase):

    def test_regular_polygons(self):
        # triangulation from a regular torus
        self.assertFalse(Triangulation.square_torus().is_non_degenerate())
        self.assertFalse(Triangulation.regular_octagon().is_non_degenerate())
        self.assertFalse(
            Triangulation.octagon_and_squares().is_non_degenerate())

    def test_arnoux_yoccoz(self):
        self.assertTrue(Triangulation.arnoux_yoccoz(3).is_non_degenerate())

    def test_shear(self):
        M = matrix([[2, 1], [1, 1]])
        sheared_torus = Triangulation.square_torus().apply_matrix(M)
        self.assertTrue(sheared_torus.is_non_degenerate())

if __name__ == "__main__":
    unittest.main(verbosity=2, failfast=True)
