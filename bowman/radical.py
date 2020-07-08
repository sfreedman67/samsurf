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


