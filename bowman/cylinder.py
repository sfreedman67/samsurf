from sage.all import *
from bowman.triangulation import Triangulation
from rational_ht_application import perp_vector_2D


class Cylinder:
    """
    Represents a Cylinder on a translation surface.

    #TODO: Only supports horizontal and vertical cylinders
    on bicuspid surfaces, for now
    """

    def __init__(self, direction, triangles, gluings):
        self.direction = direction
        self.triangles = triangles
        self.gluings = gluings

    @property
    def height(self):
        height_dir = perp_vector_2D(self.direction)
        tri = self.triangles[0]
        return tri.len_in_direction(height_dir)
        # each triangle must have an edge along cylinder direction and also
        # perpendicular direction

    @property
    def circumference(self):
        return self.width_between_tri_idx(0, len(self.triangles))

    @property
    def modulus(self):
        return self.height / self.circumference

    def width_between_tri_idx(self, start_idx, end_idx):
        """
        Returns the distance between the origin points for the triangles
        at start_idx and end_idx (with looping back allowed)
        """
        # only consider even indices, since those contain the origin point
        if start_idx % 2 == 1:
            start_idx -= 1
        if end_idx % 2 == 1:
            end_idx -= 1
        # calculate distance, skipping triangles sharing origin point
        dist_so_far = 0
        for i in range(start_idx, end_idx, 2):
            curr_idx = i % len(self.triangles)
            curr_tri = self.triangles[curr_idx]
            dist_so_far += curr_tri.len_in_direction(self.direction)
        return dist_so_far

    def get_third_constraint(self, base_tri, idx_traveled):
        """
        Finds the third constraint for the case where the periodic point starts
        in base_tri, and after veech action moves to the triangle idx_traveled
        away.
        INPUT:
        * base_tri: Triangle object
            The triangle you start from
        * idx_traveled: positive integer
            The number of triangles the periodic point travels along
        """
        start_idx = self.triangles.index(base_tri)
        end_idx = start_idx + idx_traveled
        width_travelled = self.width_between_tri_idx(start_idx, end_idx)
        other_tri = self.triangles[end_idx % len(self.triangles)]
        return other_tri.constraint_in_direction(
            self.direction, offset=-width_travelled)


    @property
    def num_twists(self, global_veech_elem):
        """
        num_twists: int
            Maximum number of times we need to mod by the width of the embedded
            triangle to return a point to its original point after applying
            the veech group element
        global_veech_elem is the parabolic element along the cylinder direction
        for the whole surface.
        """
        _, b, c, _ = *global_veech_elem[0], *global_veech_elem[1]
        # unpack matrix values
        if self.direction == vector([1, 0]):  # horizontal case
            off_diag_number = b
        elif self.direction == vector([0, 1]):  # vertical case
            off_diag_number = c
        return off_diag_number / self.modulus


"""
Old documentation:

Inputs:
* direction: 2-tuple
    Represents the direction of the cylinder in 2D space.
* height, circumference: number
    Geometric invariants of the cylinder, used in embedding
* veech_elem: 2x2 matrix
    Parabolic element corresponding to the cylinder direction
* other_dir: 2-tuple
    The direction of the other cylinders intersecting this cylinder.
* regions_dict: Dictionary, int -> (Cylinder, LinearXYPoly)
    The data about the regions this cylinder is split into by cylinders
    in the other direction. A dictionary where the keys are unique int,
    representing a specific region. The first value points to the
    cylinder in other direction this region is a part of, and the next
    value points to a LinearXYPoly representing the ratio of heights
    in the coordinates of this embedding.
"""


