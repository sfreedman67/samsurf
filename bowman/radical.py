from sage.all import *

from collections import namedtuple
from functools import lru_cache


# TODO: OK...how much speedup is *actually* coming from radical vs QQbar?
class Radical:
    def __init__(self, A, B, C):
        self.A = A
        self.B = B
        self.C = C

        self._hash = None

    def __iter__(self):
        return iter((self.A, self.B, self.C))

    @staticmethod
    @lru_cache(None)
    def sign(A, B, C):
        if B == 1:
            if C == A ** 2 and A <= 0:
                return 0
            elif C < A ** 2 and A < 0:
                return -1
            return 1
        elif B < 0:
            return -1 * Radical.sign(-A, 1, B ** 2 * C)
        else:
            return Radical.sign(A, 1, B ** 2 * C)

    @property
    def value(self):
        return QQbar(self.A) + QQbar(self.B) * QQbar(self.C).sqrt()

    def __neg__(self):
        return Radical(-self.A, -self.B, self.C)

    def __repr__(self):
        return f"{self.value}"

    def __hash__(self):
        if self._hash is None:
            self._hash = hash(self.value)
        return self._hash

    def __eq__(self, other):
        """ Checks whether A + B sqrt(C) == D + E sqrt(F)"""
        if not isinstance(other, Radical):
            return NotImplemented
        A, B, C = self
        D, E, F = other

        if Radical.sign(A - D, B, C) != Radical.sign(0, E, F):
            return False

        A1 = (A - D)**2 + B**2 * C - E**2 * F
        B1 = 2 * (A - D) * B

        return Radical.sign(A1, B1, C) == 0

    def __lt__(self, other):
        """Checks whether A + B sqrt(C) < D + E sqrt(F)"""
        if not isinstance(other, Radical):
            return NotImplemented

        A, B, C = self
        D, E, F = other

        if D != 0:
            return Radical(A - D, B, C) < Radical(0, E, F)
        elif E == 0:
            return Radical.sign(*self) < 0
        elif E == 1:
            return Radical.sign(*self) < 0 or Radical.sign(A ** 2 + B ** 2 * C - F, 2 * A * B, C) < 0
        elif E == -1:
            return Radical.sign(*self) < 0 and Radical.sign(F - A ** 2 - B ** 2 * C, -2 * A * B, C) < 0
        else:
            if E > 0:
                return self < Radical(0, 1, E ** 2 * F)
            return self < Radical(0, -1, E ** 2 * F)

    @staticmethod
    def simplify_fraction(r1, r2):
        a, b, c = r1
        d, ee, f = r2

        if c != f:
            raise ValueError("Radicals are not equal")
        elif Radical.sign(*r2) == 0:
            raise ValueError("Denominator is zero")

        if Radical.sign(d, -ee, c) == 0:
            return Radical(a / (2 * d), b / (2 * d), c)
        else:
            denom = d ** 2 - (ee ** 2) * c
            a1 = (a * d - b * ee * c) / denom
            b1 = (b * d - a * ee) / denom
            return Radical(a1, b1, c)
