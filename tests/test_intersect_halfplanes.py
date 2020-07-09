import sage.all
from sage.all import *

import unittest

import inspect

import cProfile
import pstats

from context import bowman


from bowman.intersect_halfplanes import intersect_halfplanes

from bowman.polygon import Point, Edge, plot_polygon
from bowman.radical import Radical
from bowman.halfplane import HalfPlane

from bowman.triangulation import Triangulation


class TestIntersectHalfPlanes(unittest.TestCase):

    def setUp(self):

        self.circle_neg1_pos1 = HalfPlane.from_ineq(ZZ(1), ZZ(0), ZZ(-1))
        self.circle_phi_phibar = HalfPlane.from_ineq(ZZ(-1), ZZ(1), ZZ(1))
        self.line_pos1_infty = HalfPlane.from_ineq(ZZ(0), ZZ(-1), ZZ(1))

        self.pt_neg1 = Point(ZZ(-1), ZZ(0))
        self.pt_pos1 = Point(ZZ(1), ZZ(0))
        self.pt_infty = Point(oo, ZZ(0))
        self.pt_phi = Point(Radical(QQ(1 / 2), QQ(1 / 2), QQ(5)), 0)
        self.pt_phibar = Point(Radical(QQ(1 / 2), QQ(-1 / 2), QQ(5)), 0)

        phi = Radical(QQ(1 / 2), QQ(1 / 2), QQ(5))

        self.edge_neg1_pos1 = Edge(
            self.circle_neg1_pos1, self.pt_neg1, self.pt_pos1)
        self.edge_pos1_infty = Edge(
            self.line_pos1_infty, self.pt_pos1, self.pt_infty)
        self.edge_phi_phibar = Edge(
            self.circle_phi_phibar, self.pt_phi, self.pt_phibar)

        self.ideal_pos1_neg1 = Edge(None, self.pt_pos1, self.pt_neg1)
        self.ideal_infty_pos1 = Edge(None, self.pt_infty, self.pt_pos1)
        self.ideal_infty_neg1 = Edge(None, self.pt_infty, self.pt_neg1)
        self.ideal_phibar_phi = Edge(None, self.pt_phibar, self.pt_phi)

    def test_intersecting_one_halfplane(self):
        single_halfplanes = [
            (self.circle_neg1_pos1, [
             self.edge_neg1_pos1, self.ideal_pos1_neg1]),
            (self.line_pos1_infty, [
             self.edge_pos1_infty, self.ideal_infty_pos1]),
            (self.circle_phi_phibar, [self.ideal_phibar_phi, self.edge_phi_phibar])]

        for halfplane, result in single_halfplanes:
            output = intersect_halfplanes([halfplane])
            self.assertCountEqual(output, result, msg=f"{halfplane}")

    def test_intersect_asymptotically_parallel_halfplanes(self):
        self.assertEqual(
            self.circle_neg1_pos1.intersect_boundaries(self.line_pos1_infty), self.pt_pos1)

        halfplanes = [self.circle_neg1_pos1, self.line_pos1_infty]
        result = [self.edge_neg1_pos1,
                  self.edge_pos1_infty, self.ideal_infty_neg1]

        self.assertCountEqual(intersect_halfplanes(halfplanes), result)
        # Check order doesn't matter
        self.assertCountEqual(intersect_halfplanes(halfplanes[::-1]), result)

    def test_intersect_ultraparallel_halfplanes(self):
        self.assertEqual(
            HalfPlane(1, 0, -1).intersect_boundaries(HalfPlane(1, 0, -4)), None)
        assert False, "TODO: add in"

    def test_intersect_sage_examples(self):
        examples = [(HalfPlane(1, -8, 15), HalfPlane(1, -11, 28)),
                    (HalfPlane(1, -9, 20), HalfPlane(1, -12, 35)),
                    (HalfPlane(1, 0, -1), HalfPlane(1, 0, -1))]
        answers = [Point(QQ(13 / 3), QQ(8 / 9)), Point(5, 0), None]

        for example, answer in zip(examples, answers):
            self.assertEqual(
                example[0].intersect_boundaries(example[1]), answer)

        assert False, "TODO: add in?"

    def test_intersect_AY3(self):

        X = Triangulation.arnoux_yoccoz(3)
        alpha = X.base_ring.gen()
        H = X.halfplanes()

        # X.plot_halfplanes(7)

        s0, e0 = H[0].endpoints
        s1, e1 = H[1].endpoints
        e2 = H[2].end
        s3 = H[3].start
        s4 = H[4].start
        e6 = H[6].end

        p01 = H[0].intersect_boundaries(H[1])
        p03 = H[0].intersect_boundaries(H[3])
        p34 = H[3].intersect_boundaries(H[4])
        p15 = H[1].intersect_boundaries(H[5])
        p16 = H[1].intersect_boundaries(H[6])
        p25 = H[2].intersect_boundaries(H[5])

        answers = [(),
                   (Edge(H[0], s0, e0), Edge(None, e0, s0)),
                   (Edge(H[0], s0, p01), Edge(
                       H[1], p01, e1), Edge(None, e1, s0)),
                   (Edge(H[0], s0, p01), Edge(H[1], p01, e1),
                    Edge(H[2], e1, e2), Edge(None, e2, s0)),
                   (Edge(H[3], s3, p03), Edge(H[0], p03, p01), Edge(
                       H[1], p01, e1), Edge(H[2], e1, e2), Edge(None, e2, s3)),
                   (Edge(H[1], p01, e1), Edge(H[2], e1, e2), Edge(
                       None, e2, s3), Edge(H[3], s3, p34), Edge(H[4], p34, p01)),
                   (Edge(H[1], p01, p15), Edge(H[5], p15, p25), Edge(H[2], p25, e2), Edge(
                       None, e2, s3), Edge(H[3], s3, p34), Edge(H[4], p34, p01)),
                   (Edge(None, e6, s3), Edge(H[3], s3, p34), Edge(H[4], p34, p01), Edge(H[1], p01, p16), Edge(H[6], p16, e6))]

        for idx, answer in enumerate(answers):
            output = intersect_halfplanes(H[:idx])
            self.assertCountEqual(output, answer, msg=f"\nTestCase: {idx}")

        intersect_halfplanes(H)
        # print(intersect_halfplanes(H))

        # idx = 10
        # first_k = intersect_halfplanes(H[:idx])
        # P = plot_polygon(first_k)
        # polygon_and_new_halfplane = P + H[idx].plot()
        # polygon_and_new_halfplane.save(f"first_{idx}_and_next.png")

        assert False, "Todo: Add in other partial intersections"


def run_only_one_test(name):
    suite = unittest.TestSuite()
    suite.addTest(TestIntersectHalfPlanes(name))

    runner = unittest.TextTestRunner()
    runner.run(suite)

if __name__ == "__main__":
    # unittest.main(failfast=False, verbosity=2)

    # run_only_one_test("test_intersect_AY3")

    planes_g20 = Triangulation.arnoux_yoccoz(20).halfplanes()
    cProfile.run("intersect_halfplanes(planes_g20)", "intersect.profile")
    s = pstats.Stats("intersect.profile")
    s.dump_stats("output.pstats")
    s.strip_dirs().sort_stats(pstats.SortKey.TIME).print_stats(10)
