from sage.all import *

import itertools
from itertools import combinations


class HalfPlane():

    def __init__(self, a, b, c):
        self.coefficients = (a, b, c)
        self.boundary = self._boundary()

    def __repr__(self):
        return f"HalfPlane([{self.coefficients[0]}](u^2 + v^2) + [{self.coefficients[1]}]u + {self.coefficients[2]} >= 0)"

    def __eq__(self, other):
        if isinstance(other, HalfPlane):
            return self.coefficients == other.coefficients
        return NotImplemented

    def __hash__(self):
        return hash(self.coefficients)

    def _boundary(self):
        a, b, c = self.coefficients

        discriminant = b**2 - 4 * a * c

        g = HyperbolicPlane().UHP().get_geodesic

        if bool(discriminant <= 0):
            return None
        elif a == 0:
            return g(-c / b, oo) if b < 0 else g(oo, -c / b)
        else:
            center = (-b) / (QQ(2) * a)
            radius2 = discriminant / (QQ(4) * a**2)

            radius = QQbar(radius2).sqrt()

            left = center - radius
            right = center + radius

            oriented_right = bool(a * (center**2) + b * center + c > 0)

            return g(right, left) if oriented_right else g(left, right)

    def is_interior(self, p):
        px, py = p.coordinates()
        return bool(self.coefficients[0] * (px**2 + py**2)
                    + self.coefficients[1] * px
                    + self.coefficients[2] >= 0)

    def plot(self):
        start, end = self.boundary.endpoints()
        orientation = ""

        if start != oo and end != oo:
            orientation = "orange" if start < end else "blue"

        else:
            orientation = "blue" if end == oo else "orange"

        return plot(self.boundary, axes=True, color=orientation)

    def intersection(self, h):
        g0, g1 = self.boundary, h.boundary
        if g0.is_ultra_parallel(g1):
            return None
        elif g0.is_asymptotically_parallel(g1):
            # TODO: the type of the endpts is algebraic number...
            if g0.start() in g1.endpoints():
                return (g0.start().coordinates(), 0)
            else:
                return (g0.end().coordinates(), 0)
        else:
            a, b, c = self.coefficients
            d, e, f = h.coefficients

            A = matrix([[a, b, c], [d, e, f]])
            gen_kernel = A.right_kernel().basis()[0]
            w = (1/gen_kernel[2]) * gen_kernel

            u = w[1]
            v2 = w[0] - u**2
            v = sqrt(v2)

            return (u, v)


def intersect_halfplanes(R):
    if not R:
        return []
    elif len(R) == 1:
        h0 = R[0]
        return [h0.boundary.start(), h0.boundary.end()]
    pass
