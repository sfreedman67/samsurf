#!/usr/bin/env sage

import sage.all
from sage.all import *

import collections
from collections import namedtuple


def is_nonnegative(A, B, C):
    '''returns True if A + B sqrt(C) is >= 0'''
    if bool(B == 0):
        return bool(A >= 0)

    K = (-A) / B
    # Depending on sign of b, system is either
    # b < 0 --> sqrt(C) <= K
    # b > 0 --> sqrt(C) >= K

    elif bool(B < 0):
        # sqrt(C) <= K
        return bool(K >= 0) and bool(C <= K**2)
    else:
        # sqrt(C) >= K
        return bool(K < 0) or bool(C >= K**2)


# TODO: We can't see the "oriented" attribute here, which makes some code weird
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

    def contains_coordinates(self, u, v2):
        return bool(self.a * (u**2 + v2) + self.b * u + self.c >= 0)

    def contains_start(self, other):
        raise NotImplementedError

    def intersects(self, other):
        raise NotImplementedError

    def plot(self):
        # For circles: Below Blue, Above Orange
        # For lines: Left bLue, Right oRange
        orientation = "blue" if self.oriented_positively else "orange"

        boundary = HyperbolicPlane().UHP().get_geodesic(self._start, self._end)
        return plot(boundary, axes=True, color=orientation)


class Line(HalfPlane):

    def __new__(cls, a, b, c):
        self = super(Line, cls).__new__(cls, a, b, c)
        self.endpoint_real = (-c) / b

        self.oriented_positively = bool(b < 0)

        self._start = self.endpoint_real if self.oriented_positively else oo
        self._end = oo if self.oriented_positively else self.endpoint_real

        return self

    # TODO: Modify given self is line
    def intersects(self, other):
        if isinstance(other, Line):
            return True

        '''
        Want to solve system:
        a (u^2 + v^2) + bu = -c
        d (u^2 + v^2) + eu = -f
        for "variables" (u^2 + v^2) and u
        '''

        A = matrix([[self.a, self.b], [other.a, other.b]])
        if A.determinant() == 0:
            return False

        u2_plus_v2, u = A.solve_right([-self.c, -other.c])
        v2 = u2_plus_v2 - u**2

        return (v2 >= 0)

    def _contains_line_start(self, other):
        return other._start == oo or self.contains_coordinates(other._start, 0)

    def _contains_circle_start(self, other):
        center_included = self.contains_coordinates(other.center, 0)

        d2_endpoint_to_center = (self.endpoint_real - other.center)**2
        start_closer_to_center = bool(other.radius2 <= d2_endpoint_to_center)
        endpoint_closer_to_center = bool(d2_endpoint_to_center <= other.radiu2)

        if self.oriented_positively == other.oriented_positively:
            return center_included and start_closer_to_center
        else:
            return center_included or endpoint_closer_to_center

    def contains_start(self, other):
        if isinstance(other, Line):
            return self._contains_line_start(other)

        return self._contains_circle_start(other)


class Circle(HalfPlane):

    def __new__(cls, a, b, c):
        self = super(Circle, cls).__new__(cls, a, b, c)

        self.center = (-b) / (QQ(2) * a)

        self.radius2 = (b**2 - 4 * a * c) / (QQ(4) * a**2)

        value_at_center = a * (self.center**2) + b * self.center + c
        self.oriented_positively = bool(value_at_center > 0)

        _radius = QQbar(self.radius2).sqrt()
        if self.oriented_positively:
            self._start, self._end = self.center + _radius, self.center - _radius
        else:
            self._start, self._end = self.center - _radius, self.center + _radius

        return self

    # TODO: modify given self is a circle
    def intersects(self, h):
        '''
        Want to solve system:
        a (u^2 + v^2) + bu = -c
        d (u^2 + v^2) + eu = -f
        for "variables" (u^2 + v^2) and u
        '''

        A = matrix([[self.a, self.b], [h.a, h.b]])
        if A.determinant() == 0:
            return None

        u2_plus_v2, u = A.solve_right([-self.c, -h.c])
        v2 = u2_plus_v2 - u**2

        return v2 >= 0

    def contains_start(self, other):
        pass
