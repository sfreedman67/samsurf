import sage.all
from sage.all import *

import collections
from collections import namedtuple

import functools

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
    def _is_zero(self):
        if self.B == 0 or self.C == 0:
            return self.A == 0
    
        K = (-self.A) / self.B
        return bool(K >= 0) and bool(self.C - K**2 == 0)

    @property
    def _is_negative(self):
        if self.B == 0 or self.C == 0:
            return bool(self.A < 0)
        K = (-self.A) / self.B
        if self.B > 0:
            return bool(K >= 0) and bool((self.C - K**2) < 0)
        return K <= 0 or self.C - K**2 > 0

    @property
    def value(self):
        if self._value is None:
            if self.C == 0:
                self._value = self.A
            
            self._value = self.A + self.B * AA(self.C).sqrt()
        return self._value

    def _is_valid_operand(self, other):
        return all(hasattr(other, letter) for letter in ('A', 'B', 'C'))

    def __hash__(self):
        return hash(self.value)

    def __eq__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        elif self.C == other.C == 0:
            return self.A == other.A
        return self.value  == other.value
    
    def __lt__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented

        return self.value  < other.value


