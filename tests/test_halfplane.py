import sage.all
from sage.all import *

import unittest

from context import bowman

import bowman.point_hyperbolic
from bowman.halfplane import HalfPlane
from bowman.triangulation import Triangulation


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
        self.assertFalse(
            self.unit_circle_interior.contains_point(self.i_plus_1))
        self.assertTrue(
            self.unit_circle_exterior.contains_point(self.i_plus_1))
        self.assertTrue(self.y_axis_right.contains_point(self.i_plus_1))


class TestHalfPlaneBoundary(unittest.TestCase):

    def test_endpoints_are_in_correct_order(self):
        X = Triangulation.arnoux_yoccoz(3)
        alpha = X.base_ring.gen()
        H = X.halfplanes()[6]

        self.assertEqual(H.boundary, HyperbolicPlane().UHP().get_geodesic(
            QQ(1 / 2) * alpha ** 2 + QQ(1 / 2) * alpha, QQ(-3 / 2) * alpha**2 - QQ(5 / 2) * alpha - 3))

if __name__ == "__main__":
    unittest.main(verbosity=2)
