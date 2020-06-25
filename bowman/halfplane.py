#!/usr/bin/env sage

import sage.all
from sage.all import *

import collections
from collections import namedtuple

from bowman.point_hyperbolic import Radical, Point


class HalfPlane(namedtuple('HalfPlane', ['a', 'b', 'c'])):
    __slots__ = ()

    @staticmethod
    def from_ineq(a, b, c):
        if a == 0 and b != 0:
            return Line(a, b, c)
        elif a != 0 and (b**2 - 4 * a * c) > 0:
            return Circle(a, b, c)
        else:
            raise ValueError("Coeffs determine a degenerate inequality")

    def __repr__(self):
        term_quadratic = f"[{self.a}](u^2 + v^2)+" if self.a != 0 else ""
        term_linear = f"[{self.b}]u+" if self.b != 0 else ""
        term_constant = f"[{self.c}]" if self.c != 0 else ""
        return term_quadratic + term_linear + term_constant + ">= 0"

    @property
    def oriented_positively(self):
        raise NotImplementedError

    @property
    def start(self):
        raise NotImplementedError

    @property
    def end(self):
        raise NotImplementedError

    @property
    def endpoints(self):
        return (self.start, self.end)

    def contains_point(self, point):
        if point == Point.infinity:
            return isinstance(self, Line) or (not self.oriented_positively)

        a, b, c = self
        u, v2 = point
        A, B, C = u

        # a[(A + B sqrt(C))^2 + v2] + b (A + B sqrt(C)) + c >= 0
        # --> A1 + B1 sqrt(C) >= 0

        A1 = a * (A**2 + B**2 * C + v2) + b * A + c
        B1 = a * (2 * A * B) + b * B

        return Radical(A1, B1, C).is_nonnegative

    def halfplane_intersection(self, other):
        if isinstance(self, Line) and isinstance(other, Line):
            return Point(oo, 0)

        M = matrix([[self.a, self.b], [other.a, other.b]])
        if M.determinant() == 0:
            return None

        u2_plus_v2, u = M.solve_right(vector([-self.c, -other.c]))
        v2 = u2_plus_v2 - u**2

        return None if bool(v2 < 0) else Point(Radical(u, 0, 0), v2)

    def edge_intersection(self, edge):
        return self.halfplane_intersection(edge.halfplane)

    def is_edge_exterior(self, edge):
        return not (self.contains_point(edge.start) or self.contains_point(edge.end))

    def plot(self):
        # TODO: Fix with new types for points

        # For circles: Below Blue, Above Orange
        # For lines: Left bLue, Right oRange
        orientation = "blue" if self.oriented_positively else "orange"

        # coordinate_start = oo if self.start.is_infinity else self.start.u._value
        # coordinate_end = oo if self.start.is_infinity else self.end.u._value

        value_start = self.start.u._value
        value_end = self.end.u._value

        boundary = HyperbolicPlane().UHP().get_geodesic(value_start, value_end)
        return boundary.plot(axes=True, color=orientation)


class Line(HalfPlane):

    @property
    def oriented_positively(self):
        return bool(self.b < 0)

    @property
    def start(self):
        u = -self.c / self.b if self.oriented_positively else oo
        return Point.boundary(u)

    @property
    def end(self):
        u = oo if self.oriented_positively else -self.c / self.b
        return Point.boundary(u)


class Circle(HalfPlane):

    def __new__(cls, a, b, c):
        self = super(Circle, cls).__new__(cls, a, b, c)
        self.center = Point.boundary(-b / (QQ(2) * a))
        self.radius2 = (b**2 - 4 * a * c) / (QQ(4) * a**2)

        return self

    @property
    def oriented_positively(self):
        return self.contains_point(self.center)

    @property
    def start(self):
        coord_center = self.center.u.A
        B = 1 if self.oriented_positively else -1
        r = Radical(coord_center, B, self.radius2)
        return Point.boundary(r)

    @property
    def end(self):
        coord_center = self.center.u.A
        B = -1 if self.oriented_positively else 1
        r = Radical(coord_center, B, self.radius2)
        return Point.boundary(r)


class Edge(namedtuple("Edge", ['halfplane', 'start', 'end'])):
    __slots__ = ()

    def __repr__(self):
        Ideal_descriptor = "Ideal" if self.is_ideal else ""
        return Ideal_descriptor + f"Edge from {self.start} to {self.end}"

    
    @property
    def is_ideal(self):
        return self.halfplane is not None 

    def retract_to(self, point):
        return Edge(self.halfplane, self.start, point)

    def chop_at(self, point):
        return Edge(self.halfplane, point, self.end)
