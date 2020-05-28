from sage.all import *

import itertools
from itertools import combinations

class Halfplane():

    def __init__(self, coeffs):
        self.coeffs = coeffs


def inequality_to_geodesic(a, b, c):
    g = HyperbolicPlane().UHP().get_geodesic

    if a == 0:
        if b == 0:
            raise ValueError("invalid inequality with a == b == 0")
        else:
            if b < 0:
                return g(-c / b, oo)
            else:
                return g(oo, -c / b)
    else:
        discriminant = b**2 - 4 * a * c
        if bool(discriminant <= 0):
            raise ValueError("discriminant was non-positive")
        else:
            center = (-b) / (2 * a)
            radius2 = discriminant / (4 * a**2)

            radius = AA(radius2).sqrt()

            left_root = center - radius
            right_root = center + radius

            if a * (center**2) + b * center + c > 0:
                return g(right_root, left_root)
            else:
                return g(left_root, right_root)


def is_trivial(ineq):
    return bool(ineq[1]**2 - 4 * ineq[0] * ineq[2] <= 0)


def is_interior(halfplanes, idx):
    assert False, "TODO: write"


def plot(geodesic):
    start = geodesic.start()
    end = geodesic.end()
    if not (start == oo or end == oo):
        if start < end:
            return geodesic.plot(axes=True, color="blue")
        else:
            return geodesic.plot(axes=True, color="red")
    elif start == oo:
        return geodesic.plot(axes=True, color="red")
    else:
        return geodesic.plot(axes=True, color="blue")


def intersection_pairwise(h1, h2):
    if h1.is_ultra_parallel(h2):
        return None
    elif h1.is_asymptotically_parallel(h2):
        return h1.start() if h1.start() in h2.ideal_endpoints() else h1.end()
    else:
        return h1.intersection(h2)


def intersection_halfplanes(halfplanes):
    if len(halfplanes) < 3:
        raise ValueError("Must be intersecting at least three planes")
    intersection_points = [intersection_pairwise(
        h1, h2) for h1, h2 in itertools.combinations(halfplanes, 2)]

    return intersection_points
