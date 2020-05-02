import unittest
from triangulation import Triangulation

# just want to use Sage's symbolic calculations
from sage.all import *

# input: a non-denerate Delaunay triangulation
# output: iso-dlny region

def edge_inequality(hinge):
    # want det of matrix w/ row [xi + uyi, vyi, (xi + uyi)^2 + (vyi)^2 ]
    # expand out 3rd column, remove v from second column, use multilinearity

    a = matrix([[v[0], v[1], v[1]**2] for v in hinge]).determinant()
    b = matrix([[v[0], v[1], v[0] * v[1]] for v in hinge]).determinant()
    c = matrix([[v[0], v[1], v[0]**2] for v in hinge]).determinant()

    # TODO: remove the stupid factor of 2
    # return an inequality a(u^2 + v^2) + 2bu + c >= 0
    return (a, b, c)

# Input option A
# g = (0, ctr, rad) or (1, x-coord) 
# class hyperbolicGeodesic()
# subclass vert line
# subclass circle 

# Input option B
# g = ( 1, 3) , ( 2, +inf)
# g' = ((2, 0), (4, 0))

# g = a, b . g' = c, d. What's the intersection?

# edge inequality --> 
def intersect_geodesics(g1, g2):
    # make g1 and g2 monic
    if g1[0] != 0 and g1[0] != 1:
        return intersect_geodesics((1, g1[1] / g1[0], g1[2] / g1[0]), g2)
    elif g2[0] != 0 and g2[0] != 1: 
        return intersect_geodesics(g1, (1, g2[1] / g2[0], g2[2] / g2[0]))

    # check whether g1 and g2 are the same
    if g1 == g2:
        raise ValueError("geodesics g1, g2 must be distinct")

    for g in [g1, g2]:
        # check if g1 and g2 are "valid" geodesics
        if g[0] == g[1] == 0:
            raise ValueError("geodesic cannot have a = b = 0")
        # if a != 0, make equation monic
        elif g[1]**2 - g[2] < 0:
            raise ValueError("geodesic has degenerate radius")

    a, b, c = g1[0], g1[1], g1[2]
    d, e, f = g2[0], g2[1], g2[2]

    if a == d == 0:
        # two parallel vertical lines
        return None
    # TODO: check logic
    elif a == 0 and d == 1:
        # Want the first equation to be the quadratic term
        a, d = d, a
        b, e = e, b
        c, f = f, c

    if d == 0:
        # 2nd eqn is 2eu + f = 0
        # TODO: if e = 0...
        u = -f / (2 * e)
    else:
        # solve 1st eqn for v^2, sub into 2nd eqn
        # get u^2 + (-u^2 - 2bu - c) + 2eu + f = 0
        # => u(2e - 2b) + (f - c) = 0
        # Check for e=b...
        u = (c - f) / (2 * (e - b))

    # v^2 = -u^2 -2bu - c
    # TODO: Can't take square root of negative!!
    K = -u**2 - 2 * b * u - c

    if K < 0:
        return None
    else:
        v = sqrt(K) 

    # TODO: if u and v are not both real:
    return (u, v)

def is_non_degenerate(triang):
    return all(matrix([[v[0], v[1], v[0]**2 + v[1]**2] for v in hinge]).determinant() != 0 for hinge in triang.hinges())

# input: n geodesics [(ai, bi, ci)]
# output: vertices of a hyperbolic polygon
def generate_IDR(triang):
    return None


class IsoDlnyTests(unittest.TestCase):

    def testIsNonDegenerate(self):
        # triangulation from a regular torus
        sq_torus = Triangulation.square_torus()
        self.assertFalse(is_non_degenerate(sq_torus))

        reg_oct = Triangulation.regular_octagon()
        self.assertFalse(is_non_degenerate(reg_oct))

        a_y = Triangulation.arnoux_yoccoz(3)
        self.assertTrue(is_non_degenerate(a_y))

        o_s = Triangulation.octagon_and_squares()
        self.assertFalse(is_non_degenerate(o_s))

        M = matrix([[2, 1], [1, 1]])

        sheared_torus = sq_torus.apply_matrix(M)
        self.assertTrue(is_non_degenerate(sheared_torus))

    def testEdgeInequality(self):
        h1 = (vector([2, 2]), vector([2, 4]), vector([1, 4]))
        self.assertEqual(edge_inequality(h1), (-16, -16, -4))

        # TODO: write more test cases

class IntersectingGeodesicsTests(unittest.TestCase):
    def testUltraParallelGeodesics(self):
        # two disjoint circles, two circles with one shared endpoint, etc...
        
        #two disjoint lines
        self.assertEqual(intersect_geodesics((0, 1, 1), (0, 1, -1)), None)

        #circle and disjoint line
        self.assertEqual(intersect_geodesics((1, 0, -1), (0, 1, -4)), None)
        self.assertEqual(intersect_geodesics((1, -1, -3), (0, 4, -30)), None)

    def testCircleIntersectingLine(self):
        # (a, b, c) --> (u^2 + v^2) + 2bu + c = 0 --> (u + b)^2 + v^2 = b^2 - c 
        self.assertEqual(intersect_geodesics((1, 0, -4), (0, 1/2, 0)), (0, 2))

        inf = float("+inf")
        print(1 / inf)


if __name__ == "__main__":
    unittest.main(verbosity=2)
