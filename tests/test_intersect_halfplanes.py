import unittest

from context import bowman
import bowman.triangulation as triangulation

import bowman.intersect_halfplanes
from bowman.intersect_halfplanes import intersect_halfplanes


class TestIntersectHalfPlanes(unittest.TestCase):

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

        X = triangulation.Triangulation.arnoux_yoccoz(3)
        alpha = X.base_ring.gen()
        H = X.halfplanes()
        print("start")

        s0, e0 = H[0].endpoints
        e1 = H[1].end
        e2 = H[2].end
        s3 = H[3].start

        p01 = H[0].halfplane_intersection(H[1])
        p03 = H[0].halfplane_intersection(H[3])
        p34 = H[3].halfplane_intersection(H[4])
        p15 = H[1].halfplane_intersection(H[5])
        p25 = H[2].halfplane_intersection(H[5])

        answers_intersections_partial = [[],
                                         [s0, e0],
                                         [s0, p01, e1],
                                         [s0, p01, e1, e2],
                                         [s3, p03, p01, e1, e2],
                                         [s3, p34, p01, e1, e2],
                                         [s3, p34, p01, p15, p25, e2]]

        for idx, answer in enumerate(answers_intersections_partial):
            output = intersect_halfplanes(H[:idx])
            # TODO: check whether lists are cyclic rearrangments
            self.assertEqual(output, answer, msg=f"\nTestCase: {idx}\noutput={output}\nexpected={answer}")

        print("end")

        assert False, "Todo: Add in other partial intersections"

if __name__ == "__main__":
	suite = unittest.TestSuite()
	suite.addTest(TestIntersectHalfPlanes("test_intersect_AY3"))
	runner = unittest.TextTestRunner(verbosity=1)
	runner.run(suite)
