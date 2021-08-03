from collections import namedtuple
import functools

from sage.all import *
from bowman.linear_xy_poly import LinearXYPoly


def perp_vector_2D(vec):
    """
    Returns the perpendicular vector to a 2D vector
    """
    a, b = vec
    return vector([b, -a])


def is_valid_barycentric_coordinate(a0, a1, a2):
    if a0 + a1 + a2 != 1:
        return False
    if a0 < 0 or a1 < 0 or a2 < 0:
        return False
    return True


class Triangle():
    """ A triangle with a list of marked points
    -
    :params v0, v1, v2: vectors representing edges of the triangle with property v0 + v1 + v2 = 0
    :param points_marked: an optional list containing marked points of the form
     ((a, b, c), (r, g, b)) where (a,b,c) are barycentric coords in the triangle and (r, g, b) correspond to a color for the marked point
    -
    - see a document for how the coordinates correspond to the edges
    """
    def __init__(self, v0, v1, v2, points_marked = None, lines_marked = None):
        if sum((v0, v1, v2)) != 0:
            raise ValueError("sides do not close up")
        elif sage.all.matrix([v0, -v2]).determinant() <= 0:
            raise ValueError("sides are not oriented correctly")

        self.v0 = v0
        self.v1 = v1
        self.v2 = v2

        if points_marked is None:
            self.points_marked = tuple()
        else:
            for point_marked, _ in points_marked:
                if not is_valid_barycentric_coordinate(*(point_marked)):
                    raise ValueError("Invalid barycentric coordinates.")
            self.points_marked = tuple(points_marked)

        if lines_marked is None:
            self.lines_marked = tuple()
        else:
            for start_coords, end_coords, _ in lines_marked:
                if not is_valid_barycentric_coordinate(*start_coords) or not is_valid_barycentric_coordinate(*end_coords):
                    raise ValueError("Invalid barycentric coordinates.")
            self.lines_marked = tuple(lines_marked)

    def __repr__(self):
        return f"Triangle(v0={self.v0}, v1={self.v1}, v2={self.v2}"

    def mark_point(self, coords, rgbcolor):
        """Determine if the given coordinates COORDS are valid barycentric coordinates in the Triangle self and add to points_marked if valid.
        return 0 if the cooridnates are valid, 1 otherwise
        """
        if not is_valid_barycentric_coordinate(*coords):
            raise ValueError("Invalid barycentric coordinates.")

        for i in range(len(self.points_marked)):
            point_marked, _ = self.points_marked[i]
            if coords == point_marked:
                # If the point has already been marked, just update the color.
                return Triangle(self.v0, self.v1, self.v2, self.points_marked[0:i] + ((coords, rgbcolor),) + self.points_marked[i+1:], self.lines_marked)

        # Otherwise, append the point as a new element in self.points_marked.
        return Triangle(self.v0, self.v1, self.v2, self.points_marked + ((coords, rgbcolor),), self.lines_marked)

    def mark_line(self, start_coords, end_coords, rgbcolor):
        if not is_valid_barycentric_coordinate(*start_coords) or not is_valid_barycentric_coordinate(*end_coords):
            raise ValueError("Invalid barycentric coordinates.")

        for i in range(len(self.lines_marked)):
            start_coords_marked, end_coords_marked, _ = self.lines_marked[i]
            if start_coords == start_coords_marked and end_coords == end_coords_marked:
                return Triangle(self.v0, self.v1, self.v2, self.points_marked, self.lines_marked[0:i] + ((start_coords, end_coords, rgbcolor),) + self.lines_marked[i+1:])
        try:
            tri = Triangle(self.v0, self.v1, self.v2, self.points_marked, self.lines_marked + ((start_coords, end_coords, rgbcolor),))
        except ValueError:
            print(self.v0, self.v1, self.v2)
        return Triangle(self.v0, self.v1, self.v2, self.points_marked, self.lines_marked + ((start_coords, end_coords, rgbcolor),))

    def __getitem__(self, key):
        if key == 0:
            return self.v0
        elif key == 1:
            return self.v1
        elif key == 2:
            return self.v2
        else:
            raise ValueError(f"Invalid index {key}. A triangle has only three edges.")

    def __iter__(self):
        return iter([self.v0, self.v1, self.v2])

    def reflect(self, idx):
        def reflect_vector(v, w):
            w_parallel = (v.dot_product(w) / v.dot_product(v)) * v
            w_perp = w - w_parallel
            return w - 2 * w_perp

        v_axis = self[idx]
        v_succ = self[(idx + 1) % 3]
        v_pred = self[(idx + 2) % 3]
        sides_new = {idx: -v_axis,
                     (idx + 1) % 3: reflect_vector(v_axis, -v_pred),
                     (idx + 2) % 3: reflect_vector(v_axis, -v_succ)}
        return Triangle(sides_new[0], sides_new[1], sides_new[2], self.points_marked)

    def plot(self, basepoint=sage.all.zero_vector(2)):
        # Plot the triangle.
        triangle_plot = sage.all.polygon2d(self.vertices(basepoint), fill=False).plot()

        # Plot the marked lines.
        for start_coords, end_coords, line_marked_color in self.lines_marked:
            start_cartesian = basepoint + start_coords[1]*self.v0 - start_coords[2]*self.v2
            end_cartesian = basepoint + end_coords[1]*self.v0 - end_coords[2]*self.v2
            triangle_plot = triangle_plot + sage.all.line2d([start_cartesian.numerical_approx(), end_cartesian.numerical_approx()], rgbcolor = line_marked_color).plot() 

        # Plot the marked points.
        for point_marked, point_marked_color in self.points_marked:
            point_marked_coords = basepoint + point_marked[1]*self.v0 - point_marked[2]*self.v2
            triangle_plot = triangle_plot + sage.all.point2d((float(point_marked_coords[0]), float(point_marked_coords[1])), rgbcolor = point_marked_color, size = 50).plot()
        return triangle_plot

    def vertices(self, basepoint=sage.all.zero_vector(2)):
        return [basepoint, basepoint + self.v0, basepoint - self.v2]

    def apply_matrix(self, m):
        w0 = m * self.v0
        w1 = m * self.v1
        w2 = -(w0 + w1)
        #w2 = m * self.v2
        return Triangle(w0, w1, w2, self.points_marked)

    def __hash__(self):
        return hash(tuple(coord for vertex in self.vertices() for coord in vertex))

    @property
    def area(self):
        return QQ(1/2) * abs(sage.all.matrix([self.v0, -self.v2]).determinant())

# Methods that only work for Rectilinear Triangles

    def edge_in_direction(self, direction):
        perp_dir = perp_vector_2D(direction)
        for edge in self:
            if edge.dot_product(perp_dir) == 0:
                return edge

    def len_in_direction(self, direction):
        """
        The length of the edge along direction dir
        dir is either vector([1, 0]) (horizontal) or vector([0, 1])
        (vertical)
        """
        return self.edge_in_direction(direction).norm()

    def constraint_in_direction(self, direction, offset=0):
        perp_dir = perp_vector_2D(direction)
        cyl_height = self.len_in_direction(perp_dir)
        if direction == vector([1, 0]):  # horizontal case
            return (1 / cyl_height) * LinearXYPoly([0, 1, offset])  # (y+c)/h
        elif direction == vector([0, 1]):  # vertical case
            return (1 / cyl_height) * LinearXYPoly([1, 0, offset])  # (x+c)/h

    def is_region_starter(self, direction):
        """
        returns true if the triangle is the first triangle in its region
        along the given direction
        """
        # find the vector pointing at the direction where this triangle's
        # edge should point
        if direction == vector([1, 0]):  # horizontal case
            pointing = vector([-1, 0])
        elif direction == vector([0, 1]):  # vertical case
            pointing = vector([0, 1])
        return pointing.dot_product(self.edge_in_direction(direction)) > 0
        # returns true if the edge along direction points the right way

