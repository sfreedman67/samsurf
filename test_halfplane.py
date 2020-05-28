import sage.all
from sage.all import *

import unittest
import halfplane
import triangulation
from triangulation import Triangulation
from halfplane import inequality_to_geodesic, intersection_halfplanes


class TestIntersect(unittest.TestCase):

    def setUp(self):
        self.geodesic = HyperbolicPlane().UHP().get_geodesic
        self.point = HyperbolicPlane().UHP().get_point

    def test_bad_input(self):
        for sample_input in [[], [self.geodesic(2, 3)], [self.geodesic(2, 3), self.geodesic(3, 4)]]:
            self.assertRaises(
                ValueError, halfplane.intersection_halfplanes, sample_input)

    def test_intersect_three_planes(self):
        test_cases = [([(oo, 2), (2, 3), (3, oo)], [3, oo, 2]),
                      ([(3, 6), (5, oo), (oo, 4)], [oo, sqrt(-2) + 4, sqrt(-2) + 5])]
        for (test_case, ans) in test_cases:
            self.assertEqual(intersection_halfplanes([self.geodesic(
                endpoint1, endpoint2) for endpoint1, endpoint2 in test_case]), [self.point(p) for p in ans])

    def test_input_order_doesnt_matter(self):
        assert False, "TODO: Implement me"


class TestInequalityToGeodesic(unittest.TestCase):

    def setUp(self):
        self.geodesic = HyperbolicPlane().UHP().get_geodesic

    def test_endpoints_are_in_correct_order(self):
        X = Triangulation.arnoux_yoccoz(3)
        alpha = X.gen
        E = X.edge_inequalities()[6]

        self.assertEqual(inequality_to_geodesic(*E),
                         self.geodesic(QQ(1 / 2) * alpha ** 2 + QQ(1 / 2) * alpha,
                                       QQ(-3 / 2) * alpha**2 - QQ(5 / 2) * alpha - 3))
    def test_messy_endpoints(self):
        print("starting test: messy endpoints")
        X = Triangulation.arnoux_yoccoz(3)
        alpha = X.gen
        E = X.edge_inequalities()[3]

        g = inequality_to_geodesic(*E)
        assert False, "Todo: Address the fact that taking squares are messy"

if __name__ == "__main__":
    unittest.main(verbosity=2)
