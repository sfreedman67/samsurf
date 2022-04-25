import collections
from collections import namedtuple
from functools import lru_cache

import sage.all
from sage.all import *

from samsurf import radical
from samsurf import polygon


class HalfPlane(namedtuple('HalfPlane', ['a', 'b', 'c'])):
    r""" all points u + iv such that a(u^2 + v^2) + bu + c >= 0.
    """
    __slots__ = ()

    # TODO: Why am I making these monic? Is that assumption used?
    @classmethod
    @lru_cache(None)
    def from_ineq(cls, a, b, c):
        if b ** 2 - 4 * a * c <= 0:
            raise ValueError("Coeffs are degenerate")
        elif a == 0 and b > 0:
            return Line(QQ(0), QQ(1), c / b)
        elif a == 0 and b < 0:
            return Line(QQ(0), QQ(-1), -c / b)
        elif a > 0:
            return Circle(QQ(1), b / a, c / a)
        else:
            return Circle(QQ(-1), -b / a, -c / a)

    @classmethod
    def from_two_points(cls, p0, p1):
        """
        Construct a halfplane with boundary oriented from p0 to p1
        """
        if p0 == p1:
            raise ValueError("Points are not distinct")
        elif p0 == oo or p1 == oo:
            raise ValueError("Can't (yet) build halfplane when one point is oo")
        elif p0.v2 != 0 and p1.v2 != 0:
            if p0.u.C != 0 or p1.u.C != 0:
                raise ValueError("Can (so far) only build halfplane between two interior points")
            elif p0.u == p1.u:
                return Line.from_two_points_interior(p0, p1)
            else:
                return Circle.from_two_points_interior(p0, p1)
        else:
            raise ValueError("Can (so far) only build halfplane between two interior points")

    def __repr__(self):
        term_quadratic = f"[{self.a}](u^2 + v^2)+" if self.a != 0 else ""
        term_linear = f"[{self.b}]u+" if self.b != 0 else ""
        term_constant = f"[{self.c}]" if self.c != 0 else ""
        return term_quadratic + term_linear + term_constant + ">= 0"

    @property
    def is_oriented(self):
        """
        is_oriented is True iff the shaded region is on the left of the boundary
        (with respect to its boundary orientation)
        """
        raise NotImplementedError

    @property
    def start(self):
        raise NotImplementedError

    @property
    def end(self):
        raise NotImplementedError

    @property
    def endpoints(self):
        return self.start, self.end

    @property
    def _point_inside(self):
        raise NotImplementedError

    @property
    def _point_outside(self):
        raise NotImplementedError

    @lru_cache(maxsize=None)
    def _plug_in_point(self, A, B, C, v2):
        a, b, c = self

        A1 = a * (A ** 2 + B ** 2 * C + v2) + b * A + c
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

        polygon_prev = HalfPlane.intersect_halfplanes(halfplanes[:-1])
        halfplane_curr = halfplanes[-1]

        if len(halfplanes) == 1:
            return polygon.Polygon([polygon.Edge(halfplane_curr,
                                                 halfplane_curr.start,
                                                 halfplane_curr.end),
                                    polygon.Edge(None,
                                                 halfplane_curr.end,
                                                 halfplane_curr.start)])

        elif polygon_prev is None:
            return None
        return polygon_prev.intersect_with_halfplane(halfplane_curr)

    def reorient(self):
        a, b, c = self
        return HalfPlane.from_ineq(-a, -b, -c)

    def apply_mobius(self, m):
        [[a, b], [c, d]] = m
        u, v, w = self

        u1 = (c ** 2) * w - c * d * v + (d ** 2) * u
        v1 = QQ(-2) * a * c * w + a * d * v + b * c * v - QQ(2) * b * d * u
        w1 = (a ** 2) * w - a * b * v + (b ** 2) * u

        return HalfPlane.from_ineq(u1, v1, w1)

    def plot(self):
        # For circles: Below Blue, Above Orange
        # For lines: Left bLue, Right oRange
        color_orientation = "blue" if self.is_oriented else "orange"

        value_start = oo if self.start == oo else self.start.u.value
        value_end = oo if self.end == oo else self.end.u.value

        boundary = sage.all.HyperbolicPlane().UHP().get_geodesic(value_start, value_end)
        return boundary.plot(axes=True, color=color_orientation)


class Line(HalfPlane):
    __slots__ = ()

    @property
    def is_oriented(self):
        """
        When is_oriented is true/false, line is oriented pointing north/south
        :return:
        :rtype:
        """
        return self.b == -1

    @property
    def start(self):
        if self.is_oriented:
            return polygon.Point(-self.c / self.b, QQ(0))
        return oo

    @property
    def end(self):
        if self.is_oriented:
            return oo
        return polygon.Point(-self.c / self.b, QQ(0))

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
        v2 = -(u ** 2 + b2 * u + c2)

        return None if v2 < 0 else polygon.Point(u, v2)

    def _intersect_line(self, other):
        c1 = self.c / self.b
        c2 = other.c / other.b

        if c1 == c2:
            return None
        return oo

    @classmethod
    def from_two_points_interior(cls, p0, p1):
        if p0.v2 < p1.v2:
            return Line(QQ(0), QQ(-1), p0.u.A)
        else:
            return Line(QQ(0), QQ(1), -p0.u.A)


class Circle(HalfPlane):
    __slots__ = ()

    @property
    def center(self):
        return polygon.Point(-self.b / (2 * self.a), QQ(0))

    @property
    def radius2(self):
        return (self.b ** 2 - 4 * self.a * self.c) / 4

    @property
    def is_oriented(self):
        """
        When is_oriented is true/false, circle is oriented CCW/CW
        :return:
        :rtype:
        """
        return self.a == -1

    @property
    def start(self):
        coord_center = self.center.u.A
        if self.is_oriented:
            return polygon.Point(
                radical.Radical(coord_center, QQ(1), self.radius2), QQ(0))
        else:
            return polygon.Point(
                radical.Radical(coord_center, QQ(-1), self.radius2), QQ(0))

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

    def _intersect_circle(self, other):
        if self.is_oriented:
            b1, c1 = -self.b, -self.c
        else:
            b1, c1 = self.b, self.c
        if other.is_oriented:
            b2, c2 = -other.b, -other.c
        else:
            b2, c2 = other.b, other.c

        if b1 == b2:
            return None
        else:
            u = (c2 - c1) / (b1 - b2)
            v2 = -(b1 * u + c1 + u ** 2)

            return None if v2 < 0 else polygon.Point(u, v2)

    def _intersect_line(self, other):
        return other._intersect_circle(self)

    @classmethod
    def from_two_points_interior(cls, p0, p1):
        u_0, v2_0 = p0.u.A, p0.v2
        u_1, v2_1 = p1.u.A, p1.v2
        A = sage.all.matrix([[u_0, QQ(1)], [u_1, QQ(1)]])
        Y = sage.all.vector([-u_0**2 - v2_0, -u_1**2 - v2_1])
        b, c = A.solve_right(Y)
        if p0 < p1:
            return Circle(QQ(1), b, c)
        else:
            return Circle(QQ(-1), -b, -c)
