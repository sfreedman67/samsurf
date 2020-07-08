import sage.all
from sage.all import *

import collections
from collections import namedtuple

import functools

from context import bowman
import bowman.halfplane as halfplane

@functools.total_ordering
class Radical:
    def __init__(self, A, B, C):
        if bool(C < 0):
            raise ValueError("C must be non-negative")

        self.A = A
        self.B = B
        self.C = C
        self._value = None

    def __repr__(self):
        term_rational = f"{self.A}" if bool(self.A != 0) else ""

        term_radical = "" if bool(self.B == 0) or bool(self.C == 0) else f" + {self.B}*sqrt{{{self.C}}}"

        if term_rational == "" and term_radical == "":
            return "0"
        return term_rational + term_radical

    def __iter__(self):
        return iter((self.A, self.B, self.C))

    @property
    def value(self):
        if self._value is None:
            self._value = self.A + self.B * sqrt(AA(self.C))
        return self._value

    def _is_valid_operand(self, other):
        return all(hasattr(other, letter) for letter in ('A', 'B', 'C'))

    def __hash__(self):
        return hash(self.value)

    def __eq__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented

        return self.value  == other.value
    
    def __lt__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented

        return self.value  < other.value


class Point(namedtuple('Point', ['u', 'v2'])):
    __slots__ = ()

    def __new__(cls, u, v2):
        self = super(Point, cls).__new__(cls, u, v2)

        if u != oo and not isinstance(u, Radical):
            return Point.__new__(cls, Radical(u, 0, 0), v2)

        return self

    def __repr__(self):
        return f"Point({self.u}, {self.v2})"

    @staticmethod
    def CCW(p1, p2, p3):
        if any((p1 == p2, p1 == p3, p2 == p3)):
            raise ValueError("Can only determine CCW for 3 distinct points")
        elif any(getattr(pt, "v2") != 0 for pt in (p1, p2, p3)):
            raise ValueError("Can only determine CCW for boundary points")

        if all(not pt.is_infinity for pt in (p1, p2, p3)):
            return (p1 < p2 < p3) or (p2 < p3 < p1) or (p3 < p1 < p2)

        elif p1.is_infinity:
            return p2 < p3

        elif p2.is_infinity:
            return p3 < p1

        return p1 < p2

    @property
    def is_infinity(self):
        return self.u == oo and self.v2 == 0

    @staticmethod
    def _plug_point_into_halfplane(point, plane):
        a, b, c = plane
        u, v2 = point
        A, B, C = u

        # a[(A + B sqrt(C))^2 + v2] + b (A + B sqrt(C)) + c >= 0
        # --> A1 + B1 sqrt(C) >= 0

        A1 = a * (A**2 + B**2 * C + v2) + b * A + c
        B1 = a * (2 * A * B) + b * B

        return Radical(A1, B1, C)
    # TODO: decide where this, _plug_point_into, and contains_point goes
    def is_boundary_point(self, plane):
        if self.is_infinity:
            return isinstance(plane, halfplane.Line)

        return Point._plug_point_into_halfplane(self, plane).value == 0
