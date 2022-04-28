from sage.all import *

from samsurf.polygon import Point
from samsurf.radical import Radical


def _apply_mobius_infinity(m):
    [[a, _], [c, _]] = m
    if c == 0:
        return oo
    return Point(a / c, 0)


def _apply_mobius_rational(m, pt):
    [[a, b], [c, d]] = m
    (u, _, _), v2 = pt
    if c != 0 and v2 == 0 and u == -d / c:
        return oo
    num_real = (a * u + b) * (c * u + d) + a * c * v2
    num_imag2 = v2 * (a * d - b * c)
    denom = c ** 2 * v2 + (c * u + d) ** 2

    return Point(num_real / denom, num_imag2 / denom ** 2)


def _apply_mobius_endpoint(m, pt):
    [[a, b], [c, d]] = m
    (A, B, C), _ = pt
    if c != 0 and pt.u == Radical(-d / c, 0, 0):
        return oo
    num = Radical(a * A + b, a * B, C)
    denom = Radical(c * A + d, c * B, C)
    return Point(Radical.simplify_fraction(num, denom), 0)


def apply_mobius(m, pt):
    if pt == oo:
        return _apply_mobius_infinity(m)
    elif pt.u.B == pt.u.C == 0:
        return _apply_mobius_rational(m, pt)
    return _apply_mobius_endpoint(m, pt)
