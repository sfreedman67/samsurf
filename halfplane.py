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

    def is_solution(self, p):

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
        intersections = [self.intersection_point(h) for h in halfplanes]
        intersections_valid = [point for point in intersections
                               if point is not None and in_intersection(halfplanes, point)]

        return self._order_points(_remove_duplicate_points(intersections_valid))

    def _vertices_new(self, halfplanes):
        pass



def in_intersection(halfplanes, point):
    return all(h.is_solution(point) for h in halfplanes)


def _remove_duplicate_points(points):
    coordinates = lambda p: (p.coordinates().real(), p.coordinates().imag())
    intersections_valid = list(set(map(coordinates, points)))
    return [HyperbolicPlane().UHP().get_point(x + y * QQbar(I)) for x, y in intersections_valid]


def intersect_halfplanes(halfplanes):
    vertices = []

    for i, h0 in enumerate(halfplanes):
        halfplanes_previous = halfplanes[:i]
        start, end = h0.boundary.start(), h0.boundary.end()

        vertices_new = []
        intersection_points = h0.intersection_points(halfplanes_previous)

        if in_intersection(halfplanes_previous, start) and not start in vertices:
            vertices_new.append(start)

        vertices_new.extend(
            (point for point in intersection_points if not point in vertices))

        if in_intersection(halfplanes_previous, end) and not end in vertices:
            vertices_new.append(end)

        vertices_eliminated = [
            vertex for vertex in vertices if not h0.is_solution(vertex)]

        if vertices_eliminated:
            index_min = vertices.index(vertices_eliminated[0])
            index_max = vertices.index(vertices_eliminated[-1])
            vertices = vertices[:index_min] + vertices_new + vertices[index_max + 1:]
        else:
            vertices += vertices_new

    return vertices
