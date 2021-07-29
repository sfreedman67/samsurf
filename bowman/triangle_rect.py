from sage.all import *
from linear_xy_poly import LinearXYPoly
from rational_ht_application import perp_vector_2D



class RectilinearTriangle(Triangle):
    """
    A class with methods for triangles with both a vertical
    and horizontal edge.
    Produced from bicuspid IDRs
    """

    def len_in_direction(self, direction):
        """
        The length of the edge along direction dir
        dir is either vector([1, 0]) (horizontal) or vector([0, 1])
        (vertical)
        """
        perp_dir = perp_vector_2D(direction)
        for edge in self:
            if edge.dot_product(perp_dir) == 0:
                # checks if edge is parallel to direction
                return edge.norm()

    def constraint_in_direction(self, direction, offset=0):
        cyl_height = self.len_in_direction(direction)
        if direction == vector([1, 0]):  # horizontal case
            return (1 / cyl_height) * LinearXYPoly([0, 1, offset])  # (y+c)/h
        elif direction == vector([0, 1]):  # vertical case
            return (1 / cyl_height) * LinearXYPoly([1, 0, offset])  # (x+c)/h
