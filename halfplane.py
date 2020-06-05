from sage.all import *

import itertools
from itertools import combinations

import collections
from collections import OrderedDict

import operator
from operator import itemgetter


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

    def contains_point(self, p):
        if p.coordinates() == oo:
            if self.is_circle():
                pass


        px, py = p.coordinates().real(), p.coordinates().imag()
        return bool(self.coefficients[0] * (px**2 + py**2)
                    + self.coefficients[1] * px
                    + self.coefficients[2] >= 0)

    def is_circle(self):
        return self.boundary.start() != oo and self.boundary.end() != oo

    def plot(self):
        start, end = self.boundary.endpoints()
        orientation = ""

        if self.is_circle():
            orientation = "orange" if start < end else "blue"

        else:
            orientation = "orange" if end == oo else "blue"

        return plot(self.boundary, axes=True, color=orientation)

    def intersection_point(self, h):
        # TODO: fast check for whether there is an intersection?

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
            u = QQbar(gen_ker_normalized[1])
            v2 = gen_ker_normalized[0] - u**2
            if v2 < 0:
                return None

            else:
                v = QQbar(v2).sqrt()

                return HyperbolicPlane().UHP().get_point(u + v * QQbar(I))

    def _order_points(self, points):

        if self.is_circle():
            coordinate_real = lambda p: p.coordinates().real()
            oriented_CCW = self.boundary.start() > self.boundary.end()
            return sorted(points, key=coordinate_real, reverse=oriented_CCW)
        else:
            coordinate_imag = lambda p: p.coordinates().imag()
            oriented_south = self.boundary.start() == oo
            return sorted(points, key=coordinate_imag, reverse=oriented_south)

    def intersection_points(self, halfplanes):
        points = []

        if in_halfplane_intersection(halfplanes, self.boundary.start()):
            points.append(self.boundary.start())

        valid_intersection_point = lambda point: point is not None and in_halfplane_intersection(halfplanes,
                                                                                                 point)

        intersections_boundaries = [self.intersection_point(h) for h in halfplanes
                                    if valid_intersection_point(self.intersection_point(h))]

        points.extend(intersections_boundaries)

        if in_halfplane_intersection(halfplanes, self.boundary.end()):
            points.append(self.boundary.end())

        return self._order_points(points)


def in_halfplane_intersection(halfplanes, point):
    return all(h.contains_point(point) for h in halfplanes)


def _remove_duplicate_points(points):
    coordinates = lambda p: (p.coordinates().real(), p.coordinates().imag())
    intersections_valid = list(set(map(coordinates, points)))
    return [HyperbolicPlane().UHP().get_point(x + y * QQbar(I)) for x, y in intersections_valid]


def intersect_halfplanes(halfplanes):
    vertices = []

    for idx, h0 in enumerate(halfplanes):
        vertices_new = h0.intersection_points(halfplanes[:idx])

        
        # TODO: Clunky...clean up
        indices_eliminated = [i for i, vertex in enumerate(vertices)
                              if not h0.contains_point(vertex)]

        
        if indices_eliminated:
            idx_min = indices_eliminated[0]
            idx_max = indices_eliminated[-1]
            vertices = vertices[:idx_min] + \
                vertices_new + vertices[idx_max + 1:]

        else:
            vertices += vertices_new

    return _remove_duplicate_points(vertices)
