import collections
from collections import namedtuple

import sage.all
from sage.all import *

from context import bowman
import bowman.radical
from bowman import radical
import bowman.polygon
from bowman import polygon


class HalfPlane(namedtuple('HalfPlane', ['a', 'b', 'c'])):
    __slots__ = ()

    @classmethod
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

    def _plug_in_point(self, point):
        A, B, C = point.u

        if B == 0:
            return radical.Radical(self.a * (A**2 + point.v2) + self.b * A + self.c, 0, 0)
        elif B == 1:
            return radical.Radical(self.a * (A**2 + C + point.v2) + self.b * A + self.c, 2 * self.a * A + self.b, C)

        return radical.Radical(self.a * (A**2 + C + point.v2) + self.b * A + self.c, -2 * self.a * A - self.b, C)

    def contains_point(self, point):
        if point.is_infinity:
            return isinstance(self, Line) or (not self.is_oriented)

        return not self._plug_in_point(point)._is_negative

    def contains_point_on_boundary(self, point):
        if point.is_infinity:
            return isinstance(self, Line)
        return self._plug_in_point(point)._is_zero

    def intersect_boundaries(self, other):
        if isinstance(self, Line) and isinstance(other, Line):
            return polygon.Point(oo, QQ(0))

        M = matrix([[self.a, self.b], [other.a, other.b]])
        if M.determinant() == 0:
            return None

        u2_plus_v2, u = M.solve_right(vector([-self.c, -other.c]))
        v2 = u2_plus_v2 - u**2

        return None if bool(v2 < 0) else polygon.Point(radical.Radical(u, QQ(0), QQ(0)), v2)

    def _intersect_edge_real(self, edge):
        contains_start = self.contains_point(edge.start)
        contains_end = self.contains_point(edge.end)

        if contains_start and contains_end:
            return (edge,)
        elif not (contains_start or contains_end):
            return ()

        edge_intersect_boundary = self.intersect_boundaries(edge.halfplane)

        if contains_start:
            return (polygon.Edge(edge.halfplane, edge.start, edge_intersect_boundary),)

        return (polygon.Edge(edge.halfplane,
                             edge_intersect_boundary, edge.end),)

    def _intersect_edge_ideal(self, edge):
        includes_edge_start = self.contains_point(edge.start)
        includes_edge_end = self.contains_point(edge.end)

        if includes_edge_start and includes_edge_end:
            if polygon.Point.CCW(edge.start, edge.end, self._point_outside):
                return (edge,)
            return (polygon.Edge(None, edge.start, self.start),
                    polygon.Edge(None, self.end, edge.end))
        elif includes_edge_start != includes_edge_end:
            if includes_edge_start:
                return (polygon.Edge(None, edge.start, self.start),)
            return (polygon.Edge(None, self.end, edge.end),)
        else:
            if polygon.Point.CCW(edge.start, edge.end, self._point_inside):
                return ()
            return (polygon.Edge(None, self.end, self.start),)

    def intersect_edge(self, edge):
        return self._intersect_edge_ideal(edge) if edge.is_ideal else self._intersect_edge_real(edge)

    @classmethod
    def intersect_halfplanes(cls, halfplanes):
        if not halfplanes:
            return polygon.Polygon([])

        polygon_previous = HalfPlane.intersect_halfplanes(halfplanes[:-1])
        current = halfplanes[-1]

        if polygon_previous is None:
            return None

        elif polygon_previous.edges == []:
            edges = [polygon.Edge(current, current.start, current.end),
                     polygon.Edge(None, current.end, current.start)]
            return polygon.Polygon(edges)

        return polygon_previous.intersect_with_halfplane(current)

    def plot(self):
        # For circles: Below Blue, Above Orange
        # For lines: Left bLue, Right oRange
        color_orientation = "blue" if self.is_oriented else "orange"

        value_start = oo if self.start.is_infinity else self.start.u.value
        value_end = oo if self.end.is_infinity else self.end.u.value

        boundary = HyperbolicPlane().UHP().get_geodesic(value_start, value_end)
        return boundary.plot(axes=True, color=color_orientation)


class Line(HalfPlane):
    __slots__ = ()

    @property
    def is_oriented(self):
        return bool(self.b < 0)

    @property
    def start(self):
        A = -self.c / self.b if self.is_oriented else oo
        return polygon.Point(A, QQ(0))

    @property
    def end(self):
        A = oo if self.is_oriented else -self.c / self.b
        return polygon.Point(A, QQ(0))

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


class Circle(HalfPlane):
    __slots__ = ()

    @property
    def center(self):
        return polygon.Point(-self.b / (QQ(2) * self.a), QQ(0))

    @property
    def radius2(self):
        return (self.b**2 - 4 * self.a * self.c) / QQ(4)

    @property
    def is_oriented(self):
        return self.a < 0

    @property
    def start(self):
        coord_center = self.center.u.A
        plus_or_minus = QQ(1) if self.is_oriented else QQ(-1)
        return polygon.Point(radical.Radical(coord_center, plus_or_minus, self.radius2), QQ(0))

    @property
    def end(self):
        coord_center = self.center.u.A
        plus_or_minus = QQ(-1) if self.is_oriented else QQ(1)
        return polygon.Point(radical.Radical(coord_center, plus_or_minus, self.radius2), QQ(0))

    @property
    def _point_inside(self):
        if self.is_oriented:
            return self.center
        A, B, C = self.end.u
        return polygon.Point(radical.Radical(A + 1, B, C), QQ(0))

    @property
    def _point_outside(self):
        if self.is_oriented:
            A, B, C = self.start.u
            return polygon.Point(radical.Radical(A + 1, B, C), QQ(0))
        return self.center
