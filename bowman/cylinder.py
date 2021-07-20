from sage.all import *
from bowman.triangulation import Triangulation


class Cylinder(Triangulation):
    """
    Represents a Cylinder on a translation surface.
    """

##################
# PROTOTYPE CODE #
##################

    def __init__(
        self, height, circumference, rat_ht_constraint, veech_element,
        regions={}, num_twists=1, *etc):
        """
        FOR PROTOTYPING USE ONLY
        This __init__ lets use define Cylinder objects with properties preset,
        for testing the rational height lemma calculations.
        Inputs:
        * height, circumference are geometric invariants of the cylinder
        * rat_ht_constraint is a LinearXYPoly object P. The rational height
        lemma tells us that for all periodic points (x, y) in this cylinder,
        P(x, y) is a rational number.
        * regions is a dict of int -> Cylinder object
        int represents the region's ID, and must be unique across the
        cylinder decomposition along this cylinder's direction.
        The Cylinder objects pointed to must be along the same direction.
        * veech_element is the parabolic element corresponding to this
        cylinder's direction. Can be obtained from Prop 3.4 in Wright's survey.
        * num_twists is the number of simple cut and pastes required to return
        to original cylinder shape after applying veech_element
        This is the ratio of the modulus of the cylinder with the off-diagonal
        element of veech_element
        """
        self.height = height
        self.circumference = circumference
        self.rat_ht_constraint = rat_ht_constraint
        self.regions = regions
        self.veech_element = veech_element
        self.num_twists = num_twists


###################
# INTEGRATED CODE #
###################

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

