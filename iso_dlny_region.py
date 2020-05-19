import unittest
from triangulation import Triangulation

# just want to use Sage's symbolic calculations
from sage.all import *

# input: a non-denerate Delaunay triangulation
# output: iso-dlny region


def intersect_planes(planes):
    if len(planes) < 3:
        raise ValueError("Must be intersecting at least three planes")
    return None

class TestIntersectPlanes(unittest.TestCase):
    def setUp(self):
        self.geodesic = HyperbolicPlane().UHP().get_geodesic
        A2 = matrix([[1, 2*(1 + Triangulation.regular_octagon().gen)], [0, 1]])
        self.sheared_octagon = Triangulation.regular_octagon().apply_matrix(A2)
        print(self.sheared_octagon.is_non_degenerate())

    def test_bad_input(self):
        for sample_input in [[], [self.geodesic(2, 3)], [self.geodesic(2, 3), self.geodesic(3, 4)]]:
            self.assertRaises(ValueError, intersect_planes, sample_input)

    def test_three_planes(self):
        test_cases = [([(oo, 2), (2, 3), (3, oo)], [3, oo, 2]),
                      ([(3, 6), (5, oo), (oo, 4)], [oo, sqrt(-2) + 4, sqrt(-2) + 5])]
        for (test_case, ans) in test_cases:
            self.assertEqual(intersect_planes(test_case), ans)

    def test_regular_octagon(self):
        pass # TODO add in

if __name__ == "__main__":
    unittest.main(failfast=True)
