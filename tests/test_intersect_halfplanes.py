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
        self.pt_phi = Point(Radical(QQ(1 / 2), QQ(1 / 2), QQ(5)), ZZ(0))
        self.pt_phibar = Point(Radical(QQ(1 / 2), QQ(-1 / 2), QQ(5)), ZZ(0))

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

        p68 = H[6].intersect_boundaries(H[8])
        p08 = H[0].intersect_boundaries(H[8])
        p05 = H[0].intersect_boundaries(H[5])
        p56 = H[5].intersect_boundaries(H[6])

        answer_final = [Edge(H[8], p68, p08),
                        Edge(H[0], p08, p05),
                        Edge(H[5], p05, p56),
                        Edge(H[6], p56, p68)]

        output_final = intersect_halfplanes(H)
        self.assertCountEqual(output_final, answer_final)

        hinges_from_dict = [hinge
                            for list_hinges in X.halfplanes_to_hinges.values()
                            for hinge in list_hinges]
        self.assertCountEqual(hinges_from_dict, X.hinges())

        print(f"num_distinct_halfplanes={len(X.halfplanes_to_hinges.keys())}")

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
        ineqs81 = []

        halfplanes39 = [HalfPlane.from_ineq(*ineq) for ineq in ineqs39]

        P = intersect_halfplanes(halfplanes39)

        assert False, "Todo: determine answer"


def run_only_one_test(name):
    suite = unittest.TestSuite()
    suite.addTest(TestIntersectHalfPlanes(name))

    runner = unittest.TextTestRunner()
    runner.run(suite)

if __name__ == "__main__":
    unittest.main(failfast=False, verbosity=2)

    run_only_one_test("test_regular_octagon")

    # planes = Triangulation.arnoux_yoccoz(20).halfplanes()
    # cProfile.run("intersect_halfplanes(planes)", "intersect.profile")
    # s = pstats.Stats("intersect.profile")
    # s.dump_stats("intersect.pstats")
    # s.strip_dirs().sort_stats(pstats.SortKey.TIME).print_stats(5)
