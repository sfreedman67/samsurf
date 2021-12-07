from collections import namedtuple
from functools import lru_cache

from sage.all import *

from samsurf.triangle import Triangle, is_valid_barycentric_coordinate
from samsurf.halfplane import HalfPlane


class Hinge:
    """ Two adjacent triangles where a hinge transformation could be performed.
    -
    :param tri: a Triangle
    :param id_edge: an edge vector, the shared edge identified by tri
    :param tri_opp: a Triangle, the triangle opposite tri
    :param id_edge_op: an edge vector, the shared edge identified by tri_opp
    id_edge and id_edge_op are tuples of (triangle id, edge id), not vectors

    - See document for examples of labeling tri, tri_opp and edges.
    """
    def __init__(self, tri, id_edge, tri_opp, id_edge_opp):
        self.tri = tri
        self.id_edge = id_edge
        self.tri_opp = tri_opp
        self.id_edge_opp = id_edge_opp

        edge = tri[id_edge[1]]
        edge_opp = tri_opp[id_edge_opp[1]]

        if edge != -edge_opp:
            raise ValueError("Edges either nonparallel or improperly oriented")

        v0 = tri[(id_edge[1] + 1) % 3]
        v1 = edge_opp
        v2 = -tri_opp[(id_edge_opp[1] - 1) % 3]

        self.vectors = (v0,v1,v2)

    @property
    def coordinates(self):
        return tuple(coord for vector in self.vectors for coord in vector)

    def __hash__(self):
        return hash((self.coordinates,))

    @classmethod
    def from_id_edge(cls, trin, id_edge):
        label_tri, label_edge = id_edge

        id_edge_opp = trin.gluings[id_edge]

        label_tri_opp, label_edge_opp = id_edge_opp

        tri = trin.triangles[label_tri]
        tri_opp = trin.triangles[label_tri_opp]

        return Hinge(tri, id_edge, tri_opp, id_edge_opp)

    @property
    def is_convex(self):
        v0, v1, v2 = self.vectors
        boundary = [v0, v1 - v0, v2 - v1, -v2]
        crosses = [w0x * w1y - w1x * w0y
                   for (w0x, w0y), (w1x, w1y)
                   in zip(boundary, boundary[1:] + boundary[:1])]

        all_positive = all(bool(cross > 0) for cross in crosses)
        all_negative = all(bool(cross < 0) for cross in crosses)

        return all_positive or all_negative

    @property
    def marked_cartesian_coords(self):
        v0, v1, v2 = self.vectors
        cartesian_coords_list = []

        change_of_basis_tri = column_matrix([v1, v0])
        id_edge_tri = self.id_edge[1]
        for bary_coord, color in self.tri.points_marked:
            partial_coord = (bary_coord[id_edge_tri], bary_coord[(id_edge_tri + 2) % 3])
            cartesian_coord = change_of_basis_tri * vector(partial_coord)
            cartesian_coords_list.append((cartesian_coord, color))

        change_of_basis_opp = column_matrix([v2,v1])
        id_edge_opp = self.id_edge_opp[1]
        for bary_coord, color in self.tri_opp.points_marked:
            partial_coord = (bary_coord[(id_edge_opp + 2) % 3], bary_coord[(id_edge_opp + 1) % 3])
            cartesian_coord = change_of_basis_opp * vector(partial_coord)
            cartesian_coords_list.append((cartesian_coord, color))

        return cartesian_coords_list

    def mark_cartesian_coord(self, cartesian_coord, color):
        tri = self.tri
        opp = self.tri_opp
        v0, v1, v2 = self.vectors

        change_of_basis_tri = column_matrix([v1, v0])
        id_edge_tri = self.id_edge[1]
        partial_coord = change_of_basis_tri**(-1) * vector(cartesian_coord)
        bary_coord_indexed = sorted([(id_edge_tri, partial_coord[0]),
                                     ((id_edge_tri + 1) % 3, 1 - partial_coord[0] - partial_coord[1]),
                                     ((id_edge_tri + 2) % 3, partial_coord[1])])
        bary_coord = tuple(coord for index, coord in bary_coord_indexed)
        if is_valid_barycentric_coordinate(*bary_coord):
            tri = tri.mark_point(bary_coord, color)

        change_of_basis_opp = column_matrix([v2, v1])
        id_edge_opp = self.id_edge_opp[1]
        partial_coord = change_of_basis_opp**(-1) * vector(cartesian_coord)
        bary_coord_indexed = sorted([(id_edge_opp, 1 - partial_coord[0] - partial_coord[1]),
                                     ((id_edge_opp + 1) % 3, partial_coord[1]),
                                     ((id_edge_opp + 2) % 3, partial_coord[0])])
        bary_coord = tuple(coord for index, coord in bary_coord_indexed)
        if is_valid_barycentric_coordinate(*bary_coord):
            opp = opp.mark_point(bary_coord, color)

        return Hinge(tri, self.id_edge, opp, self.id_edge_opp)

    def flip(self):
        """Performs hinge flip and maintains the locations of marked points.
        -See document for how labeling is done."""
        v0, v1, v2 = self.vectors

        # produce new side list maintaining order such that v1 is still v1
        sides_ordered = sorted([(self.id_edge[1], v0 - v2),
                                ((self.id_edge[1] + 1) % 3, v1 - v0),
                                ((self.id_edge[1] + 2) % 3, v2 - v1)])

        # produce new triangle tri from sides_ordered
        tri = Triangle(*(vector for _, vector in sides_ordered), [])

        sides_ordered = sorted([(self.id_edge_opp[1], v2 - v0),
                                ((self.id_edge_opp[1] + 1) % 3, -v2),
                                ((self.id_edge_opp[1] + 2) % 3, v0)])

        tri_opp = Triangle(*(vector for _, vector in sides_ordered), [])

        flipped_hinge = Hinge(tri, self.id_edge, tri_opp, self.id_edge_opp)

        for cartesian_coord, color in self.marked_cartesian_coords:
            flipped_hinge = flipped_hinge.mark_cartesian_coord(cartesian_coord - v0, color)

        return flipped_hinge

    @property
    def incircle_det(self):
        """(p2 is inside/on/outside oriented circle 0-P0-P1) iff (det </==/> 0) """
        return sage.all.matrix([[x, y, x ** 2 + y ** 2] for x, y in self.vectors]).determinant()

    @property
    @lru_cache(None)
    def _coefficients(self):
        (x0, y0), (x1, y1), (x2, y2) = self.vectors

        m02 = x1 * y2 - x2 * y1
        m12 = x0 * y2 - x2 * y0
        m22 = x0 * y1 - x1 * y0

        a = y0 ** 2 * m02 - y1 ** 2 * m12 + y2 ** 2 * m22
        b = 2 * (x0 * y0 * m02 - x1 * y1 * m12 + x2 * y2 * m22)
        c = x0 ** 2 * m02 - x1 ** 2 * m12 + x2 ** 2 * m22

        return a, b, c

    @property
    def halfplane(self):
        try:
            return HalfPlane.from_ineq(*self._coefficients)
        except ValueError:
            return None

    @property
    def triangle(self):
        return self.tri

    @property
    def triangle_opp(self):
        return self.tri_opp

    @property
    def ids_boundary(self):
        """return the edge IDs of the boundary of the hinge
        starting in the NE and moving Clockwise"""

        label_tri, label_edge = self.id_edge
        label_tri_opp, label_edge_opp = self.id_edge_opp

        SE = (label_tri, (label_edge + 1) % 3)
        NE = (label_tri, (label_edge + 2) % 3)
        NW = (label_tri_opp, (label_edge_opp + 1) % 3)
        SW = (label_tri_opp, (label_edge_opp + 2) % 3)

        return NE, SE, SW, NW

    def plot(self):
        v0, v1, v2 = self.vectors
        vertices_t1 = [sage.all.zero_vector(2), v0, v1]
        vertices_t2 = [sage.all.zero_vector(2), v1, v2]
        return sage.all.polygon2d(vertices_t1, fill=False).plot() + sage.all.polygon2d(vertices_t2, fill=False).plot()
