from sage.all import *
from bowman.triangulation import Triangulation


def perp_vector_2D(vec):
    """
    Returns the perpendicular vector to a 2D vector
    """
    a, b = vec
    return vector([b, -a])


class Cylinder:
    """
    Represents a Cylinder on a translation surface.

    #TODO: Only supports horizontal and vertical cylinders
    on bicuspid surfaces, for now
    """

    def __init__(self, direction, idx_dict, base_triangulation=None):
        self.direction = direction
        self.idx_dict = idx_dict  # dict: int -> triangle
        self.triangles = list(idx_dict.values())  # list of triangles
        self.base_triangulation = base_triangulation  # underlying surface

    @classmethod
    def from_indices(cls, direction, triangulation, cyl_list):
        """
        Produces a cylinder from a triangulation.
        cyl_list should be the list of indices of triangles in the cylinder
        """
        tris = [(i, triangulation.triangles[i]) for i in cyl_list]
        if tris[0][1].is_region_starter(direction):
            new_tris = tris
        else:
            new_tris = tris[1:] + tris[0:1]
        tris_dict = {idx: tri for idx, tri in new_tris}
        return Cylinder(direction, tris_dict, triangulation)

    def __iter__(self):
        return iter(self.triangles)

    def __getitem__(self, key):
        return self.triangles[key]

    @property
    def num_triangles(self):
        return len(self.triangles)

    @property
    def height(self):
        height_dir = perp_vector_2D(self.direction)
        tri = self.triangles[0]
        return tri.len_in_direction(height_dir)
        # each triangle must have an edge along cylinder direction and also
        # perpendicular direction

    @property
    def circumference(self):
        return self.width_between_tri_idx(0, self.num_triangles)

    @property
    def modulus(self):
        return self.height / self.circumference

    @property
    def other_direction(self):
        if self.direction == vector([1, 0]):
            return vector([0, 1])
        elif self.direction == vector([0, 1]):
            return vector([1, 0])
        else:
            raise ValueError(f"Cylinder direction must be horizontal"
                             f"or vertical")

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
            curr_idx = i % self.num_triangles
            curr_tri = self.triangles[curr_idx]
            dist_so_far += curr_tri.len_in_direction(self.direction)
        return dist_so_far

    def get_third_constraint(self, base_tri, idx_traveled, veech_elem):
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
        other_tri = self.triangles[end_idx % self.num_triangles]
        return other_tri.constraint_in_direction(
            self.other_direction,
            offset=-width_travelled).matrix_coords(veech_elem)

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
        answer = off_diag_number * self.modulus  # ratio with inverse of modulus
        assert int(answer) in ZZ, f"num_twists not integer with off diagonal element {off_diag_number} and modulus {self.modulus}"
        return int(answer)
