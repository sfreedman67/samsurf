import unittest

# just want to use Sage's symbolic calculations
# from sage.all import *

# input: a non-denerate Delaunay triangulation
# output: iso-dlny region


def gen_hinge(triang, edge_data1, edge_data2):
    t_lab1, t_lab2 = edge_data1[0], edge_data2[0]
    e_lab1, e_lab2 = edge_data1[1], edge_data2[1]
    tri1, tri2 = triang.triangle(t_lab1), triang.triangle(t_lab2)
    # edges are vectors
    e1, e2 = tri1.edge(e_lab1), tri2.edge(e_lab2)

    if e1 != -e2:
        raise ValueError(
            "Edges are either nonparallel or oriented incorrectly")

    v2 = e2
    v1 = tri2.edge((e_lab2 + 1) % 3)
    v3 = tri1.edge((e_lab1 - 1) % 3)

    #want p1 -> p2, p2 -> p3 to be an oriented basis
    if matrix([v_2 - v_1, v3 - v2]).determinant() < 0:
        v1, v3 = v3, v1

    return (v1, v2, v3)

def gen_edge_ineq(hinge):
    # want det of matrix w/ row [xi + uyi, vyi, (xi + uyi)^2 + (vyi)^2 ]
    # expand out 3rd column, , remove v from second column, use multilinearity

    # TODO: make cleaner
    get_coeff = lambda f: matrix([[v[0], v[1], f(v[0], v[1])]
                                  for v in hinge]).determinant()
    c = get_coeff(lambda x, y: x**2)
    b = get_coeff(operator.mul)
    a = get_coeff(lambda x, y: y**2)

    # return an inequality a(u^2 + v^2) + bu + c >= 0
    return (a, b, c)

def get_geod_intersection(g1, g2):
    # check whether g1 and g2 are the same
    if g1 == g2:
        raise ValueError("geodesics g1, g2 must be distinct")

    for g in [g1, g2]:
        # check if g1 and g2 are "valid" geodesics
        if g[0] == g[1] == 0:
            raise ValueError("geodesic cannot have a = b = 0")
        # if a != 0, make equation monic
        elif g[0] != 0:
            g[1] /= g[0]
            g[2] /= g[0]

            if g[1]**2 - g[2] <= 0:
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
    v = sqrt(-u**2 - 2 * b * u - c)

    # TODO: if u and v are not both real:
    return (u, v)

def gen_IDR(triang):
    

class IsoDlnyTests(unittest.TestCase):

    def testEdgeIneq(self):
        return None

    def testGeodesicIntersection(self):
        # two disjoint circles, two circles with one shared endpoint, etc...

if __name__ == "__main__":
    unittest.main(verbosity=2)
