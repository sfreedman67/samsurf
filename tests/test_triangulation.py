import unittest
import cProfile
import pstats

import sage.all
from sage.all import *

from context import bowman
import bowman.triangulation
from bowman import triangulation
import bowman.halfplane
from bowman import halfplane


class TestApplyMatrixToTriangle(unittest.TestCase):

    def setUp(self):
        self.torus = triangulation.Triangulation.square_torus()
        self.octagon = triangulation.Triangulation.regular_octagon()
        self.a = self.octagon.field.gen()
        self.A2 = matrix([[1, 2 * (self.a + 1)], [0, 1]])

    def test_shear_triangle0_in_torus(self):
        sq_t = triangulation.Triangulation.square_torus()
        tri0 = sq_t.triangles[0]
        M = matrix([[2, 1], [1, 1]])
        sheared_triangle = triangulation.Triangle(
            vector([3, 2]), vector([-2, -1]), vector([-1, -1]))
        self.assertEqual(tri0.apply_matrix(M), sheared_triangle)

    def test_shear_triangle5_in_octagon(self):
        tri5 = self.octagon.triangles[5]
        w1 = vector([1 / 2 * (-8 - 5 * self.a), 1 / 2 * (-2 - self.a)])
        w2 = vector([6 + 4 * self.a, 1 + self.a])
        w3 = vector([1 / 2 * (-4 - 3 * self.a), -1 / self.a])
        sheared_tri5 = triangulation.Triangle(w1, w2, w3)
        self.assertEqual(tri5.apply_matrix(self.A2), sheared_tri5)


class TestEdgeReps(unittest.TestCase):

    def test_regular_polygons(self):
        T1 = triangulation.Triangulation.square_torus()

        self.assertEqual(T1.edges, [(0, 0), (0, 1), (0, 2)])

        T2 = triangulation.Triangulation.regular_octagon()
        self.assertEqual(T2.edges, [
                         (0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (3, 0), (4, 0)])


class TestEdgeInequality(unittest.TestCase):

    def test_normal_hinge(self):
        h = triangulation.Hinge(
            (vector([2, 2]), vector([2, 4]), vector([1, 4])), (-1, -1), (-1, -1))
        self.assertEqual(h.halfplane, (-1, -2, QQ(-1 / 4)))


class TestApplyMatrixTotriangulation(unittest.TestCase):

    def test_shear_torus(self):
        sq_t = triangulation.Triangulation.square_torus()

        sheared_triangle0 = triangulation.Triangle(
            vector([3, 2]), vector([-2, -1]), vector([-1, -1]))
        sheared_triangle1 = triangulation.Triangle(
            vector([-3, -2]), vector([2, 1]), vector([1, 1]))
        sheared_triangulation = triangulation.Triangulation(
            [sheared_triangle0, sheared_triangle1], sq_t.gluings, sq_t.field)

        M = matrix([[2, 1], [1, 1]])
        self.assertEqual(sq_t.apply_matrix(M), sheared_triangulation)


class TestIsNonDegenerate(unittest.TestCase):

    def test_regular_polygons(self):
        # triangulation from a regular torus
        self.assertFalse(triangulation.Triangulation.square_torus(
        ).is_delaunay(non_degenerate=True))
        self.assertFalse(
            triangulation.Triangulation.regular_octagon().is_delaunay(non_degenerate=True))
        self.assertFalse(
            triangulation.Triangulation.octagon_and_squares().is_delaunay(non_degenerate=True))

    def test_arnoux_yoccoz(self):
        self.assertTrue(triangulation.Triangulation.arnoux_yoccoz(
            3).is_delaunay(non_degenerate=True))


class TestIsDelaunay(unittest.TestCase):

    def test_flatsurf_examples(self):
        for X in [triangulation.Triangulation.square_torus(),
                  triangulation.Triangulation.regular_octagon(),
                  triangulation.Triangulation.octagon_and_squares(),
                  triangulation.Triangulation.arnoux_yoccoz(3),
                  triangulation.Triangulation.arnoux_yoccoz(5)]:
            self.assertTrue(X.is_delaunay(non_degenerate=False))


class TestFlipHinge(unittest.TestCase):

    def test_AY3_flip_hinge(self):
        X = triangulation.Triangulation.arnoux_yoccoz(3)
        tri = X.triangles[0]
        tri_opp = X.triangles[4]

        tri_new = triangulation.Triangle(
            tri_opp[1] + tri[1], tri[2], tri_opp[0])

        tri_opp_new = triangulation.Triangle(
            tri_opp[1], tri[1], tri[2] + tri_opp[0])

        X_flipped = X.flip_hinge((0, 0))

        answer_triangles = [tri_new,
                            *X.triangles[1: 4],
                            tri_opp_new,
                            *X.triangles[5:]]

        self.assertEqual(X_flipped.triangles, answer_triangles)

        gluings_new = {(0, 0): (4, 2), (0, 1): (10, 2), (0, 2): (8, 2), (1, 0): (5, 2), (1, 1): (4, 1), (1, 2): (11, 2), (2, 0): (8, 0), (2, 1): (3, 1), (2, 2): (6, 0), (3, 0): (9, 0), (3, 1): (2, 1), (3, 2): (7, 0), (4, 0): (5, 1), (4, 1): (1, 1), (4, 2): (0, 0), (5, 0): (9, 2), (5, 1): (4, 0), (5, 2): (
            1, 0), (6, 0): (2, 2), (6, 1): (7, 1), (6, 2): (10, 0), (7, 0): (3, 2), (7, 1): (6, 1), (7, 2): (11, 0), (8, 0): (2, 0), (8, 1): (9, 1), (8, 2): (0, 2), (9, 0): (3, 0), (9, 1): (8, 1), (9, 2): (5, 0), (10, 0): (6, 2), (10, 1): (11, 1), (10, 2): (0, 1), (11, 0): (7, 2), (11, 1): (10, 1), (11, 2): (1, 2)}

        self.assertEqual(X_flipped.gluings, gluings_new)

    def test_hinge_doubly_connected(self):
        tri1 = triangulation.Triangle(vector([-1, 1]),
                                      vector([0, -2]),
                                      vector([1, 1]))

        tri2 = triangulation.Triangle(vector([1, -1]),
                                      vector([0, 2]),
                                      vector([-1, -1]))

        gluings = {(0, 0): (1, 0), (0, 1): (1, 1), (0, 2): (11, 0),
                   (1, 0): (0, 0), (1, 1): (0, 1), (1, 2): (10, 0),
                   (10, 0): (1, 2), (11, 0): (0, 2)}

        triang = triangulation.Triangulation([tri1, tri2], gluings, QQ)

        gluings_answer = {(0, 0): (10, 0), (0, 1): (1, 1), (0, 2): (1, 2),
                          (1, 0): (11, 0), (1, 1): (0, 1), (1, 2): (0, 2),
                          (10, 0): (0, 0), (11, 0): (1, 0)}

        self.assertEqual(triang.flip_hinge((0, 1)).gluings, gluings_answer)

class Test_Generate_IsoDelaunay_Complex(unittest.TestCase):
    def test_can_generate_AY3_complex(self):
        X = triangulation.Triangulation.arnoux_yoccoz(3)
        
        for num_regions in range(0, 100, 10):
            cx = X.iso_delaunay_complex(num_regions)
            fig = sum(IDR.plot() for IDR in cx)
            fig.save(f"iso_delaunay_complex_{num_regions}.png")
            
        assert False, "Todo: make a test"

if __name__ == "__main__":
    # unittest.main(verbosity=2)

    X = triangulation.Triangulation.arnoux_yoccoz(3)
    cProfile.run("X.iso_delaunay_complex(500)", "complex.profile")
    s = pstats.Stats("complex.profile")
    s.dump_stats("complex.pstats")
    s.strip_dirs().sort_stats(pstats.SortKey.CUMULATIVE).print_stats(20)
