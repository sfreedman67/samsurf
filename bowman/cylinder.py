from sage.all import *
from bowman.triangulation import Triangulation


class Cylinder:
    """
    Represents a Cylinder on a translation surface.
    """

    def __init__(self, direction, height, circumference, veech_elem,
                 other_dir, regions_dict):
        """
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
        self.direction = direction
        self.height = height
        self.circumference = circumference
        self.veech_element = veech_element
        self.other_dir = other_dir
        self.regions_dict = regions_dict

    @classmethod
    def from_triangles(self, triangles, gluings, veech_element):
        """
        Produces a Cylinder object from triangles + gluings
        """
        # compute direction, height, circumference, other_dir, regions_dict
        # from the triangles + gluings??
        # might have to pass in more info to this function
        return Cylinder(direction, height, circumference, veech_elem,
                        other_dir, regions_dict)

    @property
    def modulus(self):
        pass

    @property
    def num_twists(self):
        """
        num_twists: int
            Maximum number of times we need to mod by the width of the embedded
            triangle to return a point to its original point after applying
            the veech group element
        """
        pass

    @property
    def constraint(self):
        """
        * constraint: LinearXYPoly
            Represents the ratio of the height of a point in the cylinder with
            the cylinder's height.
        """
        pass
