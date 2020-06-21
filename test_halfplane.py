import sage.all
from sage.all import *

import unittest
import halfplane
from halfplane import HalfPlane, intersect_halfplanes
import triangulation
from triangulation import Triangulation


class TestContainsPoint(unittest.TestCase):
    def setUp(self):
        self.unit_circle_exterior = HalfPlane(1, 0, -1)
        self.unit_circle_interior = HalfPlane(-1, 0, 1)
        self.y_axis_right = HalfPlane(0, 1, 0)
        self.infty = halfplane.Point(oo, 0)
        self.i = halfplane.Point(0, 1)
        self.i_plus_1 = halfplane.Point(1, 1)

    def test_point_is_infinity(self):
        self.assertTrue(self.unit_circle_exterior.contains_point(self.infty))
        self.assertFalse(self.unit_circle_interior.contains_point(self.infty))
        self.assertTrue(self.y_axis_right.contains_point(self.infty))

    def test_point_is_finite(self):
        self.assertFalse(self.unit_circle_interior.contains_point(self.i_plus_1))
        self.assertTrue(self.unit_circle_exterior.contains_point(self.i_plus_1))
        self.assertTrue(self.y_axis_right.contains_point(self.i_plus_1))

class TestIntersectHalfPlanes(unittest.TestCase):

    def setUp(self):
        self.geodesic = HyperbolicPlane().UHP().get_geodesic
        self.point = HyperbolicPlane().UHP().get_point

    def test_intersect_asymptotically_parallel_halfplanes(self):
        unit_circle_exterior = HalfPlane(1, 0, -1)
        line_1_infty = HalfPlane(0, -1, 1)
        self.assertEqual(
            unit_circle_exterior.intersection_point(line_1_infty), (1, 0))

    def test_intersect_ultraparallel_halfplanes(self):
        self.assertEqual(
            HalfPlane(1, 0, -1).intersection_point(HalfPlane(1, 0, -4)), None)

    def test_intersect_sage_examples(self):
        examples = [(HalfPlane(1, -8, 15), HalfPlane(1, -11, 28)),
                    (HalfPlane(1, -9, 20), HalfPlane(1, -12, 35)),
                    (HalfPlane(1, 0, -1), HalfPlane(1, 0, -1))]
        answers = [(QQ(13 / 3), 2 / 3 * sqrt(2)), (5, 0), None]

        for example, answer in zip(examples, answers):
            self.assertEqual(example[0].intersection_point(example[1]), answer)

    def test_intersect_AY_initial_two_halfplanes(self):
        X = Triangulation.arnoux_yoccoz(3)
        alpha = X.base_ring.gen()
        h0, h1, *rest = X.halfplanes()

        v = QQ(1 / 2) * alpha + QQ(1 / 2)

        self.assertEqual(h0.intersection_point(h1), (0, v))

    def test_input_order_doesnt_matter(self):
        assert False, "TODO: Implement me"

    def test_intersect_AY3(self):

        X = Triangulation.arnoux_yoccoz(3)
        alpha = X.base_ring.gen()
        H = X.halfplanes()
        print("start")
        

        s0, e0 = H[0].boundary.endpoints()
        e1 = H[1].boundary.end()
        e2 = H[2].boundary.end()
        s3 = H[3].boundary.start()

        p01 = H[0].intersection_point(H[1])
        p03 = H[0].intersection_point(H[3])
        p34 = H[3].intersection_point(H[4])
        p15 = H[1].intersection_point(H[5])
        p25 = H[2].intersection_point(H[5])

        print("end")

        answers_intersections_partial = [[],
                                         [s0, e0],
                                         [s0, p01, e1],
                                         [s0, p01, e1, e2],
                                         [s3, p03, p01, e1, e2],
                                         [s3, p34, p01, e1, e2],
                                         [s3, p34, p01, p15, p25, e2]]


        for i in range(len(answers_intersections_partial)):
            self.assertEqual(intersect_halfplanes(H[:i]),
                             answers_intersections_partial[i])
        

        assert False, "Todo: Add in other partial intersections"


class TestHalfPlaneBoundary(unittest.TestCase):

    def test_endpoints_are_in_correct_order(self):
        X = Triangulation.arnoux_yoccoz(3)
        alpha = X.base_ring.gen()
        H = X.halfplanes()[6]

        self.assertEqual(H.boundary, HyperbolicPlane().UHP().get_geodesic(
            QQ(1 / 2) * alpha ** 2 + QQ(1 / 2) * alpha, QQ(-3 / 2) * alpha**2 - QQ(5 / 2) * alpha - 3))

if __name__ == "__main__":
    unittest.main(verbosity=2)
