import triangulation
from triangulation import Triangle, Triangulation
import unittest

import sage.all
from sage.all import *


class TestShearedOctagon(unittest.TestCase):

    def setUp(self):
        self.sheared_octagon = Triangulation.regular_octagon(
        ).apply_matrix(matrix(([1, 1], [0, 1])))
        self.a = Triangulation.regular_octagon().gen

    def test_can_generate_IDR(self):
        # Checks if DLNY triangulation is non-degenerate
        self.assertTrue(self.sheared_octagon.is_non_degenerate())

        # Generates the edge representatives
        self.assertEqual(self.sheared_octagon.edges(), [
                         (0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (3, 0), (4, 0)])

        # From edges, generates the hinges
        hinge1, hinge2, hinge3 = self.sheared_octagon.hinge(0, 0), self.sheared_octagon.hinge(0, 1), self.sheared_octagon.hinge(0, 2)
        self.assertEqual(hinge1, (vector(
            [self.a + 1, self.a + 1]), vector([1, 1 / 2 * self.a + 1]), vector([0, 1 / 2 * self.a])))

        self.assertEqual(hinge2, (vector(
            [-1, -1 / 2 * self.a - 1]), vector([-1, -1]), vector([-self.a - 2, -1 / 2 * self.a - 1])))

        self.assertEqual(hinge3, (vector(
            [-self.a - 2, -self.a - 1]), vector([0, -1 / 2 * self.a]), vector([1, 1])))

        # "From hinges, generates the edge inequalities"
        self.assertEqual(Triangulation.edge_inequality(hinge1), (0, self.a + 1, self.a + 1))
        self.assertEqual(Triangulation.edge_inequality(hinge2), None)
        self.assertEqual(Triangulation.edge_inequality(hinge3), (5/2 * self.a + 7/2, 0, 0))

        # From edge inequalities, generates the geodesic halfplanes
        print(self.sheared_octagon.edge_inequalities())

        # Intersects the halfplanes and forms the IDR

        assert False, "Todo: finish me"


class TestApplyMatrixToTriangle(unittest.TestCase):

    def setUp(self):
        self.torus = Triangulation.square_torus()
        self.octagon = Triangulation.regular_octagon()
        self.a = self.octagon.gen
        self.A2 = matrix([[1, 2 * (self.a + 1)], [0, 1]])

    def test_shear_triangle0_in_torus(self):
        sq_t = Triangulation.square_torus()
        tri0 = sq_t.triangles[0]
        M = matrix([[2, 1], [1, 1]])
        sheared_triangle = Triangle(
            [vector([3, 2]), vector([-2, -1]), vector([-1, -1])])
        self.assertEqual(tri0.apply_matrix(M), sheared_triangle)

    def test_shear_triangle5_in_octagon(self):
        tri5 = self.octagon.triangles[5]
        w1 = vector([1 / 2 * (-8 - 5 * self.a), 1 / 2 * (-2 - self.a)])
        w2 = vector([6 + 4 * self.a, 1 + self.a])
        w3 = vector([1 / 2 * (-4 - 3 * self.a), -1 / self.a])
        sheared_tri5 = Triangle([w1, w2, w3])
        self.assertEqual(tri5.apply_matrix(self.A2), sheared_tri5)


class TestEdgeReps(unittest.TestCase):

    def test_regular_polygons(self):
        T1 = Triangulation.square_torus()

        self.assertEqual(T1.edges(), [(0, 0), (0, 1), (0, 2)])

        T2 = Triangulation.regular_octagon()
        self.assertEqual(T2.edges(), [
                         (0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (3, 0), (4, 0)])


class TestHinge(unittest.TestCase):

    def setUp(self):
        self.octagon_and_squares = Triangulation.octagon_and_squares()

    def test_is_independent_of_edge_in_pair(self):
        for edge, opposite_edge in self.octagon_and_squares.edges(gluings=True):
            self.assertEqual(self.octagon_and_squares.hinge(
                *edge), self.octagon_and_squares.hinge(*opposite_edge))


class TestEdgeInequality(unittest.TestCase):

    def test_incorrectly_oriented_hinge(self):
        h = (vector([0, 1]), vector([1, 1]), vector([1, 0]))
        self.assertRaises(ValueError, Triangulation.edge_inequality, h)

    def test_normal_hinge(self):
        h = (vector([2, 2]), vector([2, 4]), vector([1, 4]))
        self.assertEqual(Triangulation.edge_inequality(h), (-16, -32, -4))


class TestInequaltyToGeodesic(unittest.TestCase):

    def setUp(self):
        self.geodesic = HyperbolicPlane().UHP().get_geodesic

    def test_square_torus(self):
        square_torus = Triangulation.square_torus()
        # print(square_torus.edge_inequalities())
        # [(2, 2, 0), (0, -2, 0), (0, 2, 2)]

        self.assertEqual([Triangulation.inequality_to_geodesic(ineq) for ineq in square_torus.edge_inequalities()],
                         [self.geodesic(0, oo), self.geodesic(-1, 0), self.geodesic(oo, -1)])


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
