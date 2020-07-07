import sage.all
from sage.all import *

import collections
from collections import namedtuple

from functools import total_ordering

from context import bowman
import bowman.halfplane as halfplane

@total_ordering
class Radical(namedtuple('Radical', ['A', 'B', 'C'])):
    __slots__ = ()

    def __new__(cls, A, B, C):
        self = super(Radical, cls).__new__(cls, A, B, C)

        if bool(C < 0):
            raise ValueError("C must be non-negative")

        return self

    def __repr__(self):
        term_rational = f"{self.A}" if bool(self.A != 0) else ""

        term_radical = "" if bool(self.B == 0) or bool(self.C == 0) else f" + {self.B}*sqrt{{{self.C}}}"

        if term_rational == "" and term_radical == "":
            return "0"
        return term_rational + term_radical

    def _is_valid_operand(self, other):
        return all(hasattr(other, letter) for letter in ('A', 'B', 'C'))

    def __eq__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented

        A, B, C = self
        D, E, F = other

        # A + B sqrt(C) = D + E sqrt(F)
        # (A - D) + B sqrt(C) = E sqrt(F)
        LHS = Radical(A - D, B, C)
        RHS = Radical(0, E, F)

        is_nonneg_LHS = bool(LHS.sign >= 0)
        is_nonneg_RHS = bool(RHS.sign >= 0)

        if is_nonneg_LHS != is_nonneg_RHS:
            return False

        # Square --> (A - D)^2 + B^2 C + 2(A - D)Bsqrt(C) ? E^2 F
        LHS_sq_minus_RHS_sq = Radical(
            (A - D)**2 + B**2 * C - E**2 * F, 2 * (A - D) * B, C)

        # if 0 <= x <= y , then x = y iff x^2 = y^2
        # if x <= y < 0, then x = y iff x^2 = y^2

        return LHS_sq_minus_RHS_sq.sign == 0

    def __lt__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented

        A, B, C = self
        D, E, F = other

        # A + B sqrt(C) < D + E sqrt(F)
        # (A - D) + B sqrt(C) < E sqrt(F)
        LHS = Radical(A - D, B, C)
        RHS = Radical(0, E, F)

        is_nonneg_LHS = bool(LHS.sign >= 0)
        is_nonneg_RHS = bool(RHS.sign >= 0)

        if is_nonneg_LHS and not is_nonneg_RHS:
            return False

        elif not is_nonneg_LHS and is_nonneg_RHS:
            return True

        # Square --> (A - D)^2 + B^2 C + 2(A - D)Bsqrt(C) ? E^2 F
        # [(A - D)^2 + B^2 C - E^2 F} + 2(A - D)Bsqrt(C) ? 0

        LHS_sq_minus_RHS_sq = Radical(
            (A - D)**2 + B**2 * C - E**2 * F, 2 * (A - D) * B, C)

        # if 0 <= x < y , then x < y iff x^2 < y^2
        # if x < y < 0, then x < y iff x^2 > y^2

        if is_nonneg_LHS and is_nonneg_RHS:
            return LHS_sq_minus_RHS_sq.sign < 0

        return LHS_sq_minus_RHS_sq.sign > 0

    @property
    def sign(self):
        if self.B == 0:
            return sign(self.A)

        K = -self.A / self.B

        return sign(self.B) if K < 0 else sign(self.B) * sign(self.C - K**2)


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

    def is_boundary_point(self, plane):
        if self.is_infinity:
            return isinstance(plane, halfplane.Line)

        return Point._plug_point_into_halfplane(self, plane).sign == 0

    def is_interior_point(self, plane):
        if self.is_infinity:
            return isinstance(plane, halfplane.Circle) and not (plane.is_oriented)

        return Point._plug_point_into_halfplane(self, plane).sign == 1
