import sage.all
from sage.all import *

import unittest
import halfplane
from halfplane import HalfPlane, intersect_halfplanes
import triangulation
from triangulation import Triangulation


class TestIntersectHalfPlanes(unittest.TestCase):

    def setUp(self):
        self.geodesic = HyperbolicPlane().UHP().get_geodesic
        self.point = HyperbolicPlane().UHP().get_point

    def test_intersect_asymptotically_parallel_halfplanes(self):
        unit_circle_exterior = HalfPlane(1, 0, -1)
        line_1_infty = HalfPlane(0, -1, 1)
        self.assertEqual(
            unit_circle_exterior.intersection(line_1_infty), (1, 0))

    def test_intersect_ultraparallel_halfplanes(self):
        self.assertEqual(
            HalfPlane(1, 0, -1).intersection(HalfPlane(1, 0, -4)), None)

    def test_intersect_sage_examples(self):
        examples = [(HalfPlane(1, -8, 15), HalfPlane(1, -11, 28)),
                    (HalfPlane(1, -9, 20), HalfPlane(1, -12, 35)),
                    (HalfPlane(1, 0, -1), HalfPlane(1, 0, -1))]
        answers = [(QQ(13 / 3), 2 / 3 * sqrt(2)), (5, 0), None]

        for example, answer in zip(examples, answers):
            self.assertEqual(example[0].intersection(example[1]), answer)

    def test_intersect_AY_initial_two_halfplanes(self):
        X = Triangulation.arnoux_yoccoz(3)
        alpha = X.base_ring.gen()
        h0, h1, *rest = X.halfplanes()

        v = QQ(1 / 2) * alpha + QQ(1 / 2)

        self.assertEqual(h0.intersection(h1), (0, v))

    def test_input_order_doesnt_matter(self):
        assert False, "TODO: Implement me"

    def test_intersect_AY3(self):
        print("start")

        X = Triangulation.arnoux_yoccoz(3)
        alpha = X.base_ring.gen()
        H = X.halfplanes()

        self.assertEqual(intersect_halfplanes(H[:0]), [])
        self.assertEqual(intersect_halfplanes(H[:1]), [])
        self.assertEqual(intersect_halfplanes(H[:2]), [])
        self.assertEqual(intersect_halfplanes(H[:3]), [])

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