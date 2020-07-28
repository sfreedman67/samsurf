import sage.all
from sage.all import *

import collections
from collections import namedtuple
from functools import lru_cache


class Radical(namedtuple("Radical", ["A", "B", "C"])):

    def __repr__(self):
        term_rational = f"{self.A}" if bool(self.A != 0) else ""

        term_radical = "" if bool(self.B == 0) or bool(self.C == 0) else f" + {self.B}*sqrt{{{self.C}}}"

        if term_rational == "" and term_radical == "":
            return "0"
        return term_rational + term_radical

    @property
    def _is_zero(self):
        A, B, C = self
        zero, one = A.parent().zero(), A.parent().one()

        if B == zero:
            return A == zero
        elif B == one:
            return C == A**2 and A <= zero
        else:
            return Radical(A / B, one, C)._is_zero

    @property
    def _is_negative(self):
        A, B, C = self

        zero = A.parent().zero()
        one = A.parent().one()

        if B == zero:
            return A < zero
        elif B == one:
            return A <= zero and C < A**2
        elif B == -one:
            return A < zero or A**2 < C
        else:
            if B > zero:
                return Radical(A / B, one, C)._is_negative
            return Radical(-A / B, -one, C)._is_negative

    @lru_cache(None)
    def sign(A, B, C):
        zero = A.parent().zero()
        one = A.parent().one()

        s_A = sign(A)
        s_B = sign(B)

        if s == 0:
            return s_A
        elif s == 1:
            if B == one:
                if s_A > 0:
                    return -1
                return sign(C - A**2)
            return sign(A / B, one, C)

        else:
            return -sign(-A, -B, C)

        
    @property
    def value(self):
        if self.B == 0 or self.C == 0:
            return self.A
        elif self.C.is_square():
            return self.A + self.B * self.C.sqrt()
        return QQbar(self.A) + QQbar(self.B) * QQbar(self.C).sqrt()

    def __hash__(self):
        return hash(self.value)

    def __eq__(self, other):
        if not isinstance(other, Radical):
            return NotImplemented
        A, B, C = self
        D, E, F = other

        zero, one = A.parent().zero(), A.parent().one()

        if D != zero:
            return Radical(A - D, B, C) == Radical(zero, E, F)
        elif E == zero:
            return self._is_zero
        elif E != one:
            return Radical(A / E, B / E, C) == Radical(zero, one, F)
        elif self._is_negative:
            return False
        else:
            return Radical(A**2 + B**2 * C - F, 2 * A * B, C)._is_zero

    def __lt__(self, other):
        if not isinstance(other, Radical):
            return NotImplemented

        A, B, C = self
        D, E, F = other

        zero, one = A.parent().zero(), A.parent().one()

        if D != zero:
            return Radical(A - D, B, C) < Radical(zero, E, F)
        elif E == zero:
            return self._is_negative

        elif not (E == one or E == -one):
            if E > zero:
                return Radical(A / E, B / E, C) < Radical(zero, one, F)
            return Radical(-A / E, -B / E, C) < Radical(zero, -one, F)
        elif E == one:
            return self._is_negative or Radical(A**2 + B**2 * C - F, 2 * A * B, C)._is_negative
        else:
            return not Radical(-A, -B, C)._is_negative and Radical(F - A**2 - B**2 * C, -2 * A * B, C)._is_negative
