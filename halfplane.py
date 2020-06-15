import sage.all
from sage.all import *


def ineq_is_line(a, b, c):
    return a == 0 and b != 0


def ineq_is_circle(a, b, c):
    discriminant = b**2 - 4 * a * c
    return a != 0 and discriminant > 0


def is_nondegen_ineq(a, b, c):
    return ineq_is_line(a, b, c) or ineq_is_circle(a, b, c)


class HalfPlane():

    def __init__(self, a, b, c):
        if not is_nondegen_ineq(a, b, c):
            raise ValueError("These coefficients determine a degen inequality")
        self.coefficients = (a, b, c)

    def __repr__(self):
        return f"HalfPlane([{self.coefficients[0]}](u^2 + v^2) + [{self.coefficients[1]}]u + {self.coefficients[2]} >= 0)"

    def __eq__(self, other):
        if isinstance(other, HalfPlane):
            return self.coefficients == other.coefficients
        return NotImplemented

    def __hash__(self):
        return hash(self.coefficients)

    def is_circle(self):
        return ineq_is_circle(*self.coefficients)

    def is_line(self):
        return ineq_is_line(*self.coefficients)

    def is_oriented_positively(self):
        a, b, c = self.coefficients
        if self.is_circle():
            center = (-b) / (2 * a)
            return (a * (center**2) + b * center + c) > 0
        else:
            return b < 0

    def contains_point(self, p):
        if p.coordinates() == oo:
            return self.is_line() or (self.is_circle() and self.is_oriented_positively())

        px, py = p.coordinates().real(), p.coordinates().imag()
        return bool(self.coefficients[0] * (px**2 + py**2)
                    + self.coefficients[1] * px
                    + self.coefficients[2] >= 0)

    def intersection_point(self, h):
        if self.is_line() and h.is_line():
            return oo

        '''
        Want to solve system:
        a (u^2 + v^2) + bu = -c
        d (u^2 + v^2) + eu = -f
        for "variables" (u^2 + v^2) and u
        '''
        a, b, c = self.coefficients
        d, e, f = h.coefficients

        A = matrix([[a, b], [d, e]])
        if A.determinant() == 0:
            return None

        else:
            u2_plus_v2, u = A.solve_right([-c, -f])
            v2 = u2_plus_v2 - u**2

            return None if v2 < 0 else (u, v2)

    def _order_points(self, points):
        if self.is_circle():
            coordinate_real = lambda p: p.coordinates().real()
            return sorted(points, key=coordinate_real, reverse=self.is_oriented_positively())
        else:
            coordinate_imag = lambda p: p.coordinates().imag()
            return sorted(points, key=coordinate_imag, reverse=not self.is_oriented_positively())

    def intersection_points(self, halfplanes):
        points = []

        if in_halfplane_intersection(halfplanes, self.boundary.start()):
            points.append(self.boundary.start())

        valid_intersection_point = lambda point: point is not None and in_halfplane_intersection(
            halfplanes, point)

        intersections_boundaries = [self.intersection_point(
            h) for h in halfplanes if valid_intersection_point(self.intersection_point(h))]

        points.extend(intersections_boundaries)

        if in_halfplane_intersection(halfplanes, self.boundary.end()):
            points.append(self.boundary.end())

        return self._order_points(points)

    def _boundary_geodesic(self):
        a, b, c = self.coefficients

        g = HyperbolicPlane().UHP().get_geodesic

        if self.is_line():
            return g(-c / b, oo) if self.is_oriented_positively() else g(oo, -c / b)
        else:
            center = (-b) / (QQ(2) * a)

            discriminant = b**2 - 4 * a * c
            radius2 = discriminant / (QQ(4) * a**2)
            radius = AA(radius2).sqrt()

            left = center - radius
            right = center + radius

            return g(right, left) if self.is_oriented_positively() else g(left, right)

    def plot(self):
        # For circles: Below Blue, Above Orange
        # For lines: Left bLue, Right oRange
        orientation = "blue" if self.is_oriented_positively() else "orange"

        return plot(self._boundary_geodesic(), axes=True, color=orientation)


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
