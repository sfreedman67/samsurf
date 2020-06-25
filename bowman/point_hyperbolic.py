import collections
from collections import namedtuple


class Radical(namedtuple('Radical', ['A', 'B', 'C'])):
    __slots__ = ()

    def __repr__(self):
        term_rational = f"{self.A}" if bool(self.A != 0) else ""

        term_radical = "" if bool(self.B == 0) or bool(self.C == 0) else f"{self.B}*sqrt{{{self.C}}}"

        return term_rational + term_radical

    @classmethod
    def infinity(cls):
        return Radical(oo, 0, 0)

    @property
    def is_nonnegative(self):
        '''returns True if A + B sqrt(C) is >= 0'''
        if bool(self.B == 0):
            return bool(self.A >= 0)

        K = (-self.A) / self.B
        # Depending on sign of b, system is either
        # b < 0 --> sqrt(C) <= K
        # b > 0 --> sqrt(C) >= K

        if bool(self.B < 0):
            # sqrt(C) <= K
            return bool(K >= 0) and bool(self.C <= K**2)

        # sqrt(C) >= K
        return bool(K < 0) or bool(self.C >= K**2)

    @property
    def _value(self):
        return oo if self == Radical.infinity else self.A + self.B * AA(self.C).sqrt()


class Point(namedtuple('Point', ['u', 'v2'])):
    __slots__ = ()

    def __repr__(self):
        return f"({self.u}, {self.v2})" if bool(self.v2 != 0) else f"{self.u}"

    @classmethod
    def infinity(cls):
        return Point(Radical.infinity, 0)

    @classmethod
    def boundary(cls, u):
        if not isinstance(u, Radical):
            return Point(Radical(u, 0, 0), 0)
        return Point(u, 0)

    @property
    def is_boundary(self):
        return self.v2 == 0
