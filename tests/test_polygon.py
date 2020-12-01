from unittest import TestCase

from sage.all import *

from bowman.halfplane import HalfPlane
from bowman.polygon import Edge, Point


class TestEdge(TestCase):
    def setUp(self):
        self.edge_zero_i = Edge(HalfPlane.from_ineq(0, -QQ(1), 0), Point(0, 0), Point(0, QQ(1)))
        self.edge_neg2_zero = Edge(HalfPlane.from_ineq(QQ(1), 2, 0), Point(-2, 0), Point(0, 0))

        self.edge_zero_infty = Edge(HalfPlane.from_ineq(0, -QQ(1), 0), Point(0, 0), oo)
        self.edge_infty_neg2 = Edge(HalfPlane.from_ineq(0, QQ(1), 2), oo, Point(-2, 0))

        self.edge_neg1_i = Edge(HalfPlane.from_ineq(QQ(1), 0, -QQ(1)), Point(-QQ(1), 0), Point(0, QQ(1)))
        self.edge_i_infty = Edge(HalfPlane.from_ineq(0, -QQ(1), 0), Point(0, QQ(1)), oo)
        self.edge_i_neg1 = Edge(HalfPlane.from_ineq(-QQ(1), 0, QQ(1)), Point(0, QQ(1)), Point(-QQ(1), 0))

    def test_angle_ideal_vertex(self):
        self.assertEqual(Edge.angle(self.edge_neg2_zero, self.edge_zero_i), 0)
        self.assertEqual(Edge.angle(self.edge_zero_infty, self.edge_infty_neg2), 0)

    def test_angle_interior_vertex(self):
        self.assertEqual(Edge.angle(self.edge_neg1_i, self.edge_i_infty), pi/2)
        self.assertEqual(Edge.angle(self.edge_zero_i, self.edge_i_neg1), pi/2)


