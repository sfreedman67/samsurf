# Place to store the prototype of the Rational Height Lemma

import itertools
from sage.all import *
from bowman.linear_xy_poly import LinearXYPoly
from bowman.cylinder import Cylinder


def project_polynomial(number_field, proj_matrix, poly):
    """
    Input:
        number_field is a NumberField instance,
        proj_matrix is a matrix acting on the Q-module defined by number_field
        poly is a polynomial in a polynomial ring over number_field
    Output:
        The projection of poly using proj_matrix
        i.e. proj_matrix is applied to all coefficients of poly.
    """
    coeff_dict = poly.dict()
    vector_space, vector_to_field, field_to_vector = number_field.free_module()
    for key, val in coeff_dict.items():
        vector_coeff = field_to_vector(val)
        # convert the coefficient to a vector
        projected_coeff = proj_matrix * vector_coeff
        # project the coefficient using the matrix
        coeff_dict[key] = vector_to_field(projected_coeff)
        # store the number field element corresponding to the projected vector
    return parent(poly)(coeff_dict)
    # return the polynomial corresponding to dict


def produce_segment_three_polys(list_of_poly, debug=False):
    """
    Input: list of 3 LinearXYPoly objects. Each LinearXYPoly represents linear
    P_i = a_ix + b_iy + c in Q[sqD][x, y]
    By hypothesis P_i(x, y) evaluates to a rational number
    Output: linear P such that P(x, y) = 0.
    """

    base_field = list_of_poly[0].base_field
    z_i_poly_field = PolynomialRing(base_field, 3, 'r')
    z0, z1, z2 = z_i_poly_field.gens()
    z = [z0, z1, z2]

    coeffs_matrix = [poly.get_coeffs() for poly in list_of_poly]
    # LinearXYPoly stores coeffs as a list, so this constructs matrix
    for i in range(3):
        coeffs_matrix[i][2] = z[i] - coeffs_matrix[i][2]
    Q_both = matrix(3, coeffs_matrix).determinant()
    # Combines info from all three polynomials into a single polynomial in z_i

    # separate the constraint equation into rational and irrational parts
    rat_proj_matrix = matrix(2, [[1, 0], [0, 0]])
    Q_rat = project_polynomial(base_field, rat_proj_matrix, Q_both)

    if debug:
        print("Matrix of coefficients", coeffs_matrix)
        print("Combined polynomial", Q_both)
        print("Rational polynomial", Q_rat)

    # Substitute polynomials to get linear equation
    substitution_dict = {z[i]: list_of_poly[i].get_poly() for i in range(3)}
    return Q_rat.substitute(substitution_dict)


def reduce_cylinder_constraints(cylinder, debug=False):
    """
    Finds all the segments on a cylinder where periodic points might lie.
    Inputs:
    * cylinder: Cylinder object
        Represents all the cylinders in a given direction
    Output:
        List of line segments where periodic points lie.
        Two segments are produced for each pair of regions in a given cylinder
    #TODO: This code assumes cylinder direction is horizontal, and other
    direction is vertical.
    """
    found_segments = []
    for reg_1, reg_2 in itertools.product(cylinder.regions_dict, repeat=2):
        # iterating over pairs of regions
        other_cyl_1, other_con_1 = cylinder.regions_dict[reg_1]
        other_cyl_2, other_con_2 = cylinder.regions_dict[reg_2]
        # regions_dict stores other cylinder, and embedded constraint, for a
        # given region
        if debug:
            print(
                f"Working on the possible region action "
                f"{reg_1, other_cyl_1, other_con_1} -> "
                f"{reg_2, other_cyl_1, other_con_2}")

        cylin_constraint = cylinder.constraint
        # constraint coming from cylinder cyl
        other_constraint = other_con_1
        # constraint from other cylinder that reg_1 is in
        acted_constraint = other_con_2.matrix_coords(cylinder.veech_elem)
        # constraint corresponding to new region after applying veech_elem,
        # in terms of original coordinates

        modding_offsets = [i * cylinder.circumference for i
                           in range(cylinder.num_twists + 1)]
        # possible amounts to move back by after cut and paste
        # each offset value produces its own line
        for offset in modding_offsets:
            final_constraint = acted_constraint.sub_from_x(offset)
            # constraint after moving back
            segment = produce_segment_three_polys(
                [cylin_constraint, other_constraint, final_constraint])
            found_segments.append(segment)
            #TODO: Make segment into a proper segment, rather than a polynomial
            if debug:
                print(f"constraints are {constraints}, got line {segment}")

    return found_segments







