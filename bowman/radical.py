import sage.all
from sage.all import *

import collections
from collections import namedtuple

import functools


class Radical(namedtuple("Radical", ["A", "B", "C"])):

    def __repr__(self):
        term_rational = f"{self.A}" if bool(self.A != 0) else ""

        term_radical = "" if bool(self.B == 0) or bool(self.C == 0) else f" + {self.B}*sqrt{{{self.C}}}"

        if term_rational == "" and term_radical == "":
            return "0"
        return term_rational + term_radical

    @property
    def _is_zero(self):
        if self.C == 0:
            return self.A == 0
        elif self.B == 0:
            return self.A == 0
        elif self.C * self.B**2 - self.A**2 != 0:
            return False
        elif self.B > 0:
            return self.A <= 0
        return self.A >= 0


    @property
    def _is_negative(self):
        A, B, C = self
        if B == 0:
            return A < 0
        elif not (B == QQ(1) or B == QQ(-1)):
            if B > 0:
                return Radical(A / B, QQ(1), C)._is_negative
            return Radical(-A / B, QQ(-1), C)._is_negative
        elif B == QQ(1):
            return A <= 0 and C < A**2
        else:                
            return A < 0 or A**2 < C

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
        A, B, C = self
        D, E, F = other

        if D != QQ(0):
            return Radical(A - D, B, C) == Radical(QQ(0), E, F)
        elif E == QQ(0):
            return self._is_zero
        elif E != QQ(1):
            return Radical(A / E, B / E, C) == Radical(QQ(0), QQ(1), F)
        elif self._is_negative:
            return False
        else:
            return Radical(A**2 + B**2 * C - F, 2 * A * B, C)._is_zero

    def __lt__(self, other):
        A, B, C = self
        D, E, F = other

        if E == 0:
            return Radical(A - D, B, C)._is_negative

        G = (A - D) / E
        H = B // E

        LHS = Radical(G, H, C)

        if E > 0:
            LHS_sq_minus_F = Radical(G**2 + H**2 * C - F, 2 * G * H, C)
            return LHS._is_negative or LHS_sq_minus_F._is_negative
        else:
            F_minus_LHS_sq = Radical(F - G**2 - H**2 * C, -2 * G * H, C)
            return not LHS._is_negative and F_minus_LHS_sq._is_negative
