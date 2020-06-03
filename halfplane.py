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

            radius = AA(radius2).sqrt()

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
            orientation = "orange" if end == oo else "blue"

        return plot(self.boundary, axes=True, color=orientation)

    def intersect_boundaries(self, h):
        if self.boundary.is_ultraparallel(h.boundary):
            return None
        elif self.boundary.is_asympototically_parallel(h.boundary):
            if self.boundary.start() in h.boundary.endpoints():
                return self.boundary.start()
            else:
                return self.boundary.end()
        else:
            return self.boundary.intersection(h.boundary)[0]

    def intersection(self, h):
        ''' return *the* point of intersection of boundary of self and h'''
        if self == h:
            return None

        a, b, c = self.coefficients
        d, e, f = h.coefficients

        A = matrix([[a, b, c], [d, e, f]])
        gen_ker = A.right_kernel().basis()[0]

        if gen_ker[2] == 0:
            return None
        else:
            gen_ker_normalized = (1 / gen_ker[2]) * gen_ker
            u = AA(gen_ker_normalized[1])
            v2 = gen_ker_normalized[0] - u**2
            if v2 < 0:
                return None
            
            else:
                v = AA(v2).sqrt()

                return (u, v)


def intersect_halfplanes(R):
    if not R:
        return []
    
    return [h0.boundary.intersection(h1.boundary) for h0, h1 in itertools.combinations(R, 2)]
