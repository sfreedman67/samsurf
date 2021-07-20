from sage.all import *
from bowman.triangulation import Triangulation

class Cylinder(Triangulation):
    """
    Represents a Cylinder on a translation surface.
    """

    # TODO: maybe updated __init__ with circumference and height calculated?

    def get_direction(self):
        """
        What direction is the cylinder in?
        Is this something that must be handed in through __init__?
        """
        pass

    def get_circumference(self):
        """
        Compute the circumference of the cylinder, from the triangles + gluings
        IDEA: use straight line flow in cylinder direction to compute the 
        distance traveled, 
        """
        pass

    def get_height(self):
        """
        Compute the height of the cylinder, from the triangles + gluings
        IDEA: Get the height to edge along cyl direction in one triangle
        """
        pass

    def get_modulus(self):
        pass

    def get_veech_element(self, parent_triangulation):
        """
        parent_triangulation represents the Translation Surface cylinder is on
        Returns the minimum parabolic element corresponding to this cylinder's
        direction.
        """
        pass

