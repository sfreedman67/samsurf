import unittest
import cProfile
import pstats

import sage.all
from sage.all import *

from context import bowman
import bowman.polygon
from bowman import polygon
import bowman.radical
from bowman import radical
import bowman.halfplane
from bowman import halfplane
import bowman.triangulation
from bowman import triangulation


class TestIntersectHalfPlanes(unittest.TestCase):

    def setUp(self):

        self.circle_neg1_pos1 = halfplane.HalfPlane.from_ineq(
            ZZ(1), ZZ(0), ZZ(-1))
        self.circle_phi_phibar = halfplane.HalfPlane.from_ineq(
            ZZ(-1), ZZ(1), ZZ(1))
        self.line_pos1_infty = halfplane.HalfPlane.from_ineq(
            ZZ(0), ZZ(-1), ZZ(1))

        self.pt_neg1 = polygon.Point(ZZ(-1), ZZ(0))
        self.pt_pos1 = polygon.Point(ZZ(1), ZZ(0))
        self.pt_infty = polygon.Point(oo, ZZ(0))
        self.pt_phi = polygon.Point(radical.Radical(
            QQ(1 / 2), QQ(1 / 2), QQ(5)), ZZ(0))
        self.pt_phibar = polygon.Point(
            radical.Radical(QQ(1 / 2), QQ(-1 / 2), QQ(5)), ZZ(0))

        phi = radical.Radical(QQ(1 / 2), QQ(1 / 2), QQ(5))

        self.edge_neg1_pos1 = polygon.Edge(
            self.circle_neg1_pos1, self.pt_neg1, self.pt_pos1)
        self.edge_pos1_infty = polygon.Edge(
            self.line_pos1_infty, self.pt_pos1, self.pt_infty)
        self.edge_phi_phibar = polygon.Edge(
            self.circle_phi_phibar, self.pt_phi, self.pt_phibar)

        self.ideal_pos1_neg1 = polygon.Edge(None, self.pt_pos1, self.pt_neg1)
        self.ideal_infty_pos1 = polygon.Edge(None, self.pt_infty, self.pt_pos1)
        self.ideal_infty_neg1 = polygon.Edge(None, self.pt_infty, self.pt_neg1)
        self.ideal_phibar_phi = polygon.Edge(None, self.pt_phibar, self.pt_phi)

    def test_intersecting_one_halfplane(self):
        single_halfplanes = [(self.circle_neg1_pos1, [self.edge_neg1_pos1, self.ideal_pos1_neg1]),
                             (self.line_pos1_infty, [self.edge_pos1_infty, self.ideal_infty_pos1]),
                             (self.circle_phi_phibar, [self.ideal_phibar_phi, self.edge_phi_phibar])]

        for plane, result in single_halfplanes:
            output = halfplane.HalfPlane.intersect_halfplanes([plane])
            self.assertCountEqual(output.edges, result, msg=f"{halfplane}")

    def test_intersect_asymptotically_parallel_halfplanes(self):
        self.assertEqual(
            self.circle_neg1_pos1.intersect_boundaries(self.line_pos1_infty), self.pt_pos1)

        halfplanes = [self.circle_neg1_pos1, self.line_pos1_infty]
        result = [self.edge_neg1_pos1,
                  self.edge_pos1_infty, self.ideal_infty_neg1]

        self.assertCountEqual(
            halfplane.HalfPlane.intersect_halfplanes(halfplanes).edges, result)
        # Check order doesn't matter
        self.assertCountEqual(halfplane.HalfPlane.intersect_halfplanes(
            halfplanes[::-1]).edges, result)

    def test_intersect_ultraparallel_halfplanes(self):
        h0 = halfplane.HalfPlane.from_ineq(QQ(1), 0, QQ(-1))
        h1 = halfplane.HalfPlane.from_ineq(QQ(1), 0, QQ(-4))
        h1bar = halfplane.HalfPlane.from_ineq(QQ(-1), 0, QQ(4))

        s0, e0 = h0.endpoints
        s1, e1 = h1.endpoints

        self.assertCountEqual(
            halfplane.HalfPlane.intersect_halfplanes([h0, h1]).edges,
            polygon.Polygon([polygon.Edge(None, e1, s1),
                             polygon.Edge(h1, s1, e1)]).edges
        )

        self.assertCountEqual(
            halfplane.HalfPlane.intersect_halfplanes([h1, h0]).edges,
            polygon.Polygon([polygon.Edge(None, e1, s1),
                             polygon.Edge(h1, s1, e1)]).edges
        )

        self.assertCountEqual(
            halfplane.HalfPlane.intersect_halfplanes([h0, h1bar]).edges,
            polygon.Polygon([polygon.Edge(None, s1, s0), polygon.Edge(h0, s0, e0),
                             polygon.Edge(None, e0, e1), polygon.Edge(h1bar, e1, s1)]).edges
        )

    def test_intersect_sage_examples(self):
        h00, h01 = (halfplane.HalfPlane.from_ineq(QQ(1), QQ(-8), QQ(15)),
                    halfplane.HalfPlane.from_ineq(QQ(1), QQ(-11), QQ(28)))
        p00_01 = h00.intersect_boundaries(h01)
        s00 = h00.start
        e01 = h01.end

        self.assertCountEqual(halfplane.HalfPlane.intersect_halfplanes([h00, h01]).edges,
                              [polygon.Edge(h00, s00, p00_01),
                               polygon.Edge(h01, p00_01, e01),
                               polygon.Edge(None, e01, s00)
                               ])

        h10, h11 = [halfplane.HalfPlane.from_ineq(QQ(1), QQ(-9), QQ(20)),
                    halfplane.HalfPlane.from_ineq(QQ(1), QQ(-12), QQ(35))]
        s10, e10 = h10.endpoints
        s11, e11 = h11.endpoints

        p01 = h10.intersect_boundaries(h11)

        self.assertEqual(p01, e10)
        self.assertEqual(p01, s11)

        self.assertCountEqual(halfplane.HalfPlane.intersect_halfplanes([h10, h11]).edges,
                              [polygon.Edge(h10, s10, e10),
                               polygon.Edge(h11, s11, e11),
                               polygon.Edge(None, e11, s10)
                               ])

        h20 = halfplane.HalfPlane.from_ineq(1, 0, -1)
        s20, e20 = h20.endpoints

        self.assertCountEqual(halfplane.HalfPlane.intersect_halfplanes([h20, h20]).edges,
                              [polygon.Edge(h20, s20, e20), polygon.Edge(None, e20, s20)])

    def test_intersect_AY3(self):
        X = triangulation.Triangulation.arnoux_yoccoz(3)
        alpha = X.base_ring.gen()

        H = X.halfplanes

        p68 = H[6].intersect_boundaries(H[8])
        p08 = H[0].intersect_boundaries(H[8])
        p05 = H[0].intersect_boundaries(H[5])
        p56 = H[5].intersect_boundaries(H[6])

        answer_final = [polygon.Edge(H[8], p68, p08),
                        polygon.Edge(H[0], p08, p05),
                        polygon.Edge(H[5], p05, p56),
                        polygon.Edge(H[6], p56, p68)]

        output_final = halfplane.HalfPlane.intersect_halfplanes(H)
        self.assertCountEqual(output_final.edges, answer_final)

        hinges_from_dict = [hinge
                            for list_hinges in X.halfplanes_to_hinges.values()
                            for hinge in list_hinges]
        self.assertCountEqual(hinges_from_dict, X.hinges)

    def test_regular_octagon(self):
        K = QuadraticField(2)
        a = K.gen()

        ineqs39 = [(QQ(1 / 2) * a - 1, a - 2, QQ(-1 / 2) * a + 1),
                   (2 * a - 2, 4 * a - 6, 0),
                   (-a, -2 * a, a),
                   (0, 2 * a, 0),
                   (2 * a + 2, -2, 0),
                   (4 * a + 4, -4, 0),
                   (2 * a - 2, -4 * a + 6, 0)]

        halfplanes39 = list(set(halfplane.HalfPlane.from_ineq(*ineq)
                                for ineq in ineqs39))

        p01 = halfplanes39[0].intersect_boundaries(halfplanes39[1])
        p03 = halfplanes39[0].intersect_boundaries(halfplanes39[3])
        p13 = halfplanes39[1].intersect_boundaries(halfplanes39[3])

        self.assertCountEqual(
            halfplane.HalfPlane.intersect_halfplanes(halfplanes39).edges,
            polygon.Polygon([polygon.Edge(halfplanes39[0], p01, p03),
                             polygon.Edge(halfplanes39[3], p03, p13),
                             polygon.Edge(halfplanes39[1], p13, p01)]).edges
        )


def run_only_one_test(name):
    suite = unittest.TestSuite()
    suite.addTest(TestIntersectHalfPlanes(name))

    runner = unittest.TextTestRunner()
    runner.run(suite)

if __name__ == "__main__":
    unittest.main(failfast=False, verbosity=2)

    # run_only_one_test("test_regular_octagon")

    # planes = triangulation.Triangulation.arnoux_yoccoz(20).halfplanes()
    # cProfile.run("halfplane.HalfPlane.intersect_halfplanes(planes)", "intersect.profile")
    # s = pstats.Stats("intersect.profile")
    # s.dump_stats("intersect.pstats")
    # s.strip_dirs().sort_stats(pstats.SortKey.TIME).print_stats(5)
