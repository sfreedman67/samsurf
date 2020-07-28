import collections
from collections import namedtuple
from functools import lru_cache

import sage.all
from sage.all import *

from context import bowman
import bowman.radical
from bowman import radical
import bowman.polygon
from bowman import polygon
import bowman.intersect_convex_polygons
from bowman import intersect_convex_polygons


class HalfPlane(namedtuple('HalfPlane', ['a', 'b', 'c'])):
    __slots__ = ()

    @classmethod
    @lru_cache(None)
    def from_ineq(cls, a, b, c):
        if b**2 - 4 * a * c <= 0:
            raise ValueError("Coeffs are degenerate")
        elif a == 0 and b > 0:
            return Line(0, 1, c / b)
        elif a == 0 and b < 0:
            return Line(0, -1, -c / b)
        elif a > 0:
            return Circle(1, b / a, c / a)
        else:
            return Circle(-1, -b / a, -c / a)

    def __repr__(self):
        term_quadratic = f"[{self.a}](u^2 + v^2)+" if self.a != 0 else ""
        term_linear = f"[{self.b}]u+" if self.b != 0 else ""
        term_constant = f"[{self.c}]" if self.c != 0 else ""
        return term_quadratic + term_linear + term_constant + ">= 0"

    @property
    def is_oriented(self):
        raise NotImplementedError

    @property
    def start(self):
        raise NotImplementedError

    @property
    def end(self):
        raise NotImplementedError

    @property
    def endpoints(self):
        return (self.start, self.end)

    @property
    def _point_inside(self):
        raise NotImplementedError

    @property
    def _point_outside(self):
        raise NotImplementedError

    @lru_cache(maxsize=None)
    def _plug_in_point(self, A, B, C, v2):
        a, b, c = self

        A1 = a * (A**2 + B**2 * C + v2) + b * A + c
        B1 = a * 2 * A * B + b * B
        C1 = C

        return radical.Radical(A1, B1, C1)

    def contains_point(self, point):
        if point == oo:
            return isinstance(self, Line) or (not self.is_oriented)

        output = self._plug_in_point(*point.u, point.v2)
        return radical.Radical.sign(*output) >= 0

    def contains_point_on_boundary(self, point):
        if point == oo:
            return isinstance(self, Line)
        output = self._plug_in_point(*point.u, point.v2)
        return radical.Radical.sign(*output) == 0

    def intersect_boundaries(self, other):
        if isinstance(other, Circle):
            return self._intersect_circle(other)

        return self._intersect_line(other)

    def _intersect_edge_real(self, edge, contains_start, contains_end):
        if contains_start and contains_end:
            return (edge,)
        elif not (contains_start or contains_end):
            return ()

        edge_intersect_boundary = self.intersect_boundaries(edge.halfplane)

        if contains_start:
            return (polygon.Edge(edge.halfplane, edge.start, edge_intersect_boundary),)

        return (polygon.Edge(edge.halfplane,
                             edge_intersect_boundary, edge.end),)

    def _intersect_edge_ideal(self, edge, contains_start, contains_end):
        if contains_start and contains_end:
            if polygon.Point.CCW(edge.start, edge.end, self._point_outside):
                return (edge,)
            return (polygon.Edge(None, edge.start, self.start),
                    polygon.Edge(None, self.end, edge.end))
        elif contains_start and not contains_end:
            return (polygon.Edge(None, edge.start, self.start),)
        elif not contains_start and contains_end:
            return (polygon.Edge(None, self.end, edge.end),)
        else:
            if polygon.Point.CCW(edge.start, edge.end, self._point_inside):
                return ()
            return (polygon.Edge(None, self.end, self.start),)

    def intersect_edge(self, edge, contains_start, contains_end):
        if edge.is_ideal:
            return self._intersect_edge_ideal(edge, contains_start, contains_end)

        return self._intersect_edge_real(edge, contains_start, contains_end)

    @staticmethod
    def intersect_halfplanes(halfplanes):
        if not halfplanes:
            return polygon.Polygon([])

        polygon_current = polygon.Polygon([])

        for idx in range(len(halfplanes)):
            current = halfplanes[idx]

            if polygon_current is None:
                return None

            elif idx == 0:
                edges = [polygon.Edge(current, current.start, current.end),
                         polygon.Edge(None, current.end, current.start)]
                polygon_current = polygon.Polygon(edges)

            else:
                polygon_current = polygon_current.intersect_with_halfplane(
                    current)

        return polygon_current

    def reorient(self):
        a, b, c = self
        return HalfPlane.from_ineq(-a, -b, -c)

    def plot(self):
        # For circles: Below Blue, Above Orange
        # For lines: Left bLue, Right oRange
        color_orientation = "blue" if self.is_oriented else "orange"

        value_start = oo if self.start == oo else self.start.u.value
        value_end = oo if self.end == oo else self.end.u.value

        boundary = HyperbolicPlane().UHP().get_geodesic(value_start, value_end)
        return boundary.plot(axes=True, color=color_orientation)


class Line(HalfPlane):
    __slots__ = ()

    @property
    def is_oriented(self):
        return self.b == -1

    @property
    def start(self):
        if self.is_oriented:
            return polygon.Point(-self.c / self.b, 0)
        return oo

    @property
    def end(self):
        if self.is_oriented:
            return oo
        return polygon.Point(-self.c / self.b, 0)

    @property
    def endpoint_real(self):
        return self.start if self.is_oriented else self.end

    @property
    def _point_inside(self):
        A, *_ = self.endpoint_real.u
        if self.is_oriented:
            return polygon.Point(A - 1, 0)
        return polygon.Point(A + 1, 0)

    @property
    def _point_outside(self):
        A, *_ = self.endpoint_real.u
        if self.is_oriented:
            return polygon.Point(A + 1, 0)
        return polygon.Point(A - 1, 0)

    def _intersect_circle(self, other):
        c1 = self.c / self.b
        b2, c2 = other.b / other.a, other.c / other.a

        u = -c1
        v2 = -(u**2 + b2 * u + c2)

        return None if v2 < 0 else polygon.Point(u, v2)

    def _intersect_line(self, other):
        c1 = self.c / self.b
        c2 = other.c / other.b

        if c1 == c2:
            return None
        return oo


class Circle(HalfPlane):
    __slots__ = ()

    @property
    def center(self):
        return polygon.Point(-self.b / (2 * self.a), 0)

    @property
    def radius2(self):
        return (self.b**2 - 4 * self.a * self.c) / 4

    @property
    def is_oriented(self):
        return self.a == -1

    @property
    def start(self):
        coord_center = self.center.u.A
        if self.is_oriented:
            return polygon.Point(
                radical.Radical(coord_center, 1, self.radius2), 0)
        else:
            return polygon.Point(
                radical.Radical(coord_center, -1, self.radius2), 0)

    @property
    def end(self):
        coord_center = self.center.u.A
        plus_or_minus = -1 if self.is_oriented else 1
        return polygon.Point(radical.Radical(coord_center, plus_or_minus, self.radius2), 0)

    @property
    def _point_inside(self):
        if self.is_oriented:
            return self.center
        A, B, C = self.end.u
        return polygon.Point(radical.Radical(A + 1, B, C), 0)

    @property
    def _point_outside(self):
        if self.is_oriented:
            A, B, C = self.start.u
            return polygon.Point(radical.Radical(A + 1, B, C), 0)
        return self.center

    def _intersect_circle(self, other):
        b1, c1 = self.b / self.a, self.c / self.a
        b2, c2 = other.b / other.a, other.c / other.a

        if b1 == b2:
            return None
        else:
            u = (c2 - c1) / (b1 - b2)
            v2 = -(b1 * u + c1 + u**2)

            return None if v2 < 0 else polygon.Point(u, v2)

    def _intersect_line(self, other):
        return other._intersect_circle(self)
