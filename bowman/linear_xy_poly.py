# Defines a class LinearXYPoly, which has all the properties we want
# for linear poly ax+by+c

from sage.all import *
from sage.rings.number_field.number_field import is_QuadraticField


def common_quadratic_field(list_of_nums):
    # returns QQ if all the numbers in list_of_nums is rational
    # otherwise returns the quadratic field they're all in
    output = QQ
    found_quadratic = False
    for num in list_of_nums:
        if num in QQ:
            continue
        else:
            num_field = parent(num)
            if is_QuadraticField(num_field):
                if not found_quadratic:
                    output = num_field  # update when first quadratic field found
                    found_quadratic = True
                elif num_field.discriminant() == output.discriminant():
                    continue
                else:
                    raise ValueError(f"Input {num} in different field from {output}")
            else:
                raise ValueError(f"Inputs must be in quadratic field, which {num} is not.")
    return output


class LinearXYPoly:
    """
    A linear polynomial in x, y
    Can be initialized from either list of coefficients or polynomial object
    #TODO: the base field must be QQ or a Quadratic extension
    """

    def __init__(self, coeffs_list, base_field=None):
        """
        Input must be the list of coefficients [a,b,c]
        """

        if isinstance(coeffs_list, list):  # if input is list
            assert len(coeffs_list) == 3, "Input must have 3 coefficients"
        else:
            raise ValueError("Input is not list of coefficients or polynomial")

        if base_field is None:
            self.base_field = common_quadratic_field(coeffs_list)
        else:
            self.base_field = base_field
        self.coeffs = tuple([self.base_field(x) for x in coeffs_list])
        self.poly_ring = PolynomialRing(self.base_field, 2, 'xy')
        x, y = self.poly_ring.gens()
        self.x, self.y = x, y

    @classmethod
    def from_polynomial(cls, polynomial):

        self.poly_ring = parent(polynomial)
        x, y = self.poly_ring.gens()
        self.x, self.y = x, y
        self.base_field = self.poly_ring.base_ring()

        assert polynomial.degree(x) in [0, 1], "Polynomial not linear in x"
        assert polynomial.degree(y) in [0, 1], "Polynomial not linear in y"
        coeffs_list = [
            polynomial.coefficient([1, 0]),
            polynomial.coefficient([0, 1]),
            polynomial.coefficient([0, 0])]
        # outputs coefficient of x, y, 1
        self.coeffs = tuple([self.base_field(x) for x in coeffs_list])


    def __str__(self):
        a, b, c = self.coeffs
        return f"LinearXYPoly Object representing {a}*x + {b}*y + {c}."

    def __repr__(self):
        return f"LinearXYPoly({self.get_coeffs()})"

    def can_deal_with_num(self, num):
        # tries to deal with num in self's base field, throws error if it
        # can't be dealth with
        if num in self.base_field:
            return
        elif self.base_field == QQ and is_QuadraticField(parent(num)):
            self.base_field = parent(num)
            return
        else:
            raise ValueError(f"{num} is not in polynomial's base field.")
            return

    def __sub__(self, num):
        """
        num is a number. Returns the polynomial corresponding to self,
        with other subtracted from constant term.
        """
        can_deal_with_num(self, num)
        a, b, c = self.coeffs
        return LinearXYPoly([a, b, c - num], base_field=self.base_field)

    def __rmul__(self, num):
        """
        num is a number. Returns the polynomial with num multiplied to
        the coefficients
        """
        self.can_deal_with_num(num)
        a, b, c = self.coeffs
        new_coeffs = [num * a, num * b, num * c]
        return LinearXYPoly(new_coeffs, base_field=self.base_field)

    def get_poly(self):
        """
        Outputs the polynomial associated with self.
        """
        a, b, c = self.coeffs
        return a*self.x + b*self.y + c

    def get_coeffs(self):
        return list(self.coeffs)

    def matrix_coords(self, mat):
        """
        Outputs the polynomial obtained by first applying matrix mat to x, y,
        and then applying self on the new coordinates.
        Treats x, y as a basis and apply the 2x2 matrix to the polynomial.
        Example:
        * poly = ax + by + c
        * matrix = [[p, q], [r, s]]
        * output = (ap + br)x + (aq + bs)y + c
        """
        a, b, c = self.coeffs
        p, q, r, s = *mat[0], *mat[1]  # unpack matrix [[p, q], [r, s]]
        [self.can_deal_with_num(num) for num in (p, q, r, s)]
        return LinearXYPoly([a*p + b*r, a*q + b*s, c],
                            base_field=self.base_field)



"""
TODO: Other things that would be nice to have:
Left multiplication by numbers -> multiply all coefficients by the number
"""
