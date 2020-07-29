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

    @staticmethod
    @lru_cache(None)
    def sign(A, B, C):
        if B == 1:
            if C == A**2 and A <= 0:
                return 0
            elif C < A**2 and A < 0:
                return -1
            return 1
        elif B < 0:
            return -1 * Radical.sign(-A, 1, B**2 * C)
        else:
            return Radical.sign(A, 1, B**2 * C)

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
        ''' Checks whether A + B sqrt(C) == D + E sqrt(F)'''
        if not isinstance(other, Radical):
            return NotImplemented
        A, B, C = self
        D, E, F = other

        if D != 0:
            return Radical(A - D, B, C) == Radical(0, E, F)
        elif E == 1:
            return Radical.sign(*self) >= 0 and Radical.sign(A**2 + B**2 * C - F, 2 * A * B, C) == 0
        elif E < 0:
            return Radical(-A, -B, C) == Radical(0, 1, E**2 * F)
        else:
            return self == Radical(0, 1, E**2 * F)

    def __lt__(self, other):
        '''Checks whether A + B sqrt(C) < D + E sqrt(F)'''
        if not isinstance(other, Radical):
            return NotImplemented

        A, B, C = self
        D, E, F = other

        if D != 0:
            return Radical(A - D, B, C) < Radical(0, E, F)
        elif E == 0:
            return Radical.sign(*self) < 0
        elif E == 1:
            return Radical.sign(*self) < 0 or Radical.sign(A**2 + B**2 * C - F, 2 * A * B, C) < 0
        elif E == -1:
            return Radical.sign(*self) < 0 and Radical.sign(F - A**2 - B**2 * C, -2 * A * B, C) < 0
        else:
            if E > 0:
                return self < Radical(0, 1, E**2 * F)
            return self < Radical(0, -1, E**2 * F)
