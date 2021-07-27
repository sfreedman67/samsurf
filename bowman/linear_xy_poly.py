# Defines a class LinearXYPoly, which has all the properties we want
# for linear poly ax+by+c

from sage.all import *


class LinearXYPoly:
    """
    A linear polynomial in x, y
    Can be initialized from either list of coefficients or polynomial object
    """

    def __init__(self, input_val, base_field=QQ):
        """
        Input can be the list of coefficients [a,b,c]
        Or a linear polynomial object in k[x,y]
        """
        self.base_field = base_field
        self.poly_ring = PolynomialRing(base_field, 2, 'xy')
        x, y = self.poly_ring.gens()
        self.x, self.y = x, y

        if isinstance(input_val, list):  # if input is list
            assert len(input_val) == 3, "Input must have 3 coefficients"
            self.coeffs = tuple(input_val)

        elif isinstance(input_val, type(x + y + 1)):  # if input is polynomial
            assert input_val.degree(x) in [0, 1], "Polynomial not linear in x"
            assert input_val.degree(y) in [0, 1], "Polynomial not linear in y"
            coeffs_list = [
                input_val.coefficient([1, 0]),
                input_val.coefficient([0, 1]),
                input_val.coefficient([0, 0])]
            # outputs coefficient of x, y, 1
            self.coeffs = tuple([base_field(x) for x in coeffs_list])

        else:
            raise ValueError("Input is not list of coefficients or polynomial")

    def __str__(self):
        a, b, c = self.coeffs
        return f"LinearXYPoly Object representing {a}*x + {b}*y + {c}."

    def __repr__(self):
        return f"LinearXYPoly({self.get_coeffs()})"

    def __sub__(self, num):
        """
        num is a number. Returns the polynomial corresponding to self,
        with other subtracted from constant term.
        """
        a, b, c = self.coeffs
        return LinearXYPoly([a, b, c - num], base_field=self.base_field)

    def get_poly(self):
        """
        Outputs the polynomial associated with self.
        """
        a, b, c = self.coeffs
        return a*self.x + b*self.y + c

    def get_coeffs(self):
        return list(self.coeffs)

    def apply_matrix(self, matrix):
        """
        Treats x, y as a basis and apply the 2x2 matrix to the polynomial.
        Output is LinearXYPoly object with matrix applied.
        Example:
        * poly = ax + by + c
        * matrix = [[p, q], [r, s]]
        * output = (ap + br)x + (aq + bs)y + c
        """
        a, b, c = self.coeffs
        p, q, r, s = *matrix[0], *matrix[1]  # unpack matrix [[p, q], [r, s]]
        return LinearXYPoly([a*p + b*r, a*q + b*s, c],
                            base_field=self.base_field)



"""
TODO: Other things that would be nice to have:
Left multiplication by numbers -> multiply all coefficients by the number
Left addition by addition -> add to constant term?
"""
