# Place to store the prototype of the Rational Height Lemma

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


def produce_segment_three_polys(list_of_poly):
    """
    Input: list of 3 LinearXYPoly objects. Each LinearXYPoly represents linear
    P_i = a_ix + b_iy + c in Q[sqD][x, y]
    By hypothesis P_i(x, y) evaluates to a rational number
    Output: linear P such that P(x, y) = 0.
    """

    # manually adding the variables 
    #TODO: fix for general number of polynomials

    base_field = list_of_poly[0].base_field
    r_i_poly_field = PolynomialRing(base_field, 3, 'r')
    r0, r1, r2 = r_i_poly_field.gens()
    r = [r0, r1, r2]

    coeffs_matrix = [poly.coeffs for poly in list_of_poly]
    # LinearXYPoly stores coeffs as a list, so this constructs matrix
    for i in range(3):
        coeffs_matrix[i][2] -= r[i]
    Q_both = matrix(3, coeffs_matrix).determinant()
    # Combines info from all three polynomials into a single polynomial in r_i

    # separate the constraint equation into rational and irrational parts

    rat_proj_matrix = matrix(2, [[1, 0], [0, 0]])
    Q_rat = project_polynomial(
        base_field, rat_proj_matrix, Q_both)
    assert (
        Q_rat != R(0),
        "Degenerate case encountered, projection polynomial identically zero."
    )  # Check for degenerate case. If encountered, try other directions

    irr_proj_matrix = matrix(2, [[0, 0], [0, 1]])
    Q_irr = project_polynomial(
        base_field, irr_proj_matrix, big_constraint)
    # produce the irrational equation,
    # not used since it's the polynomial  -(rat_constraint)

    # Substitute polynomials to get linear equation
    substitution_dict = {r[i]: list_of_poly[i].get_poly() for i in range(3)}
    return Q_rat.substitute(substitution_dict)

# constraints are tied to each cylinder
# also the output constraint is different from the input constraint
# code these in!!


def reduce_cylinder_constraints(cyl_decomposition, debug=False):
    """
    Inputs:
        cyl_decomposition is a list of Cylinder objects.
        refer to cylinder.py, PROTOTYPE CODE, for more info.
    Output:
        List of line segments where periodic points lie.
        Two segments are produced for each pair of regions in a given cylinder
    Assumes the veech group element being applied is horizontal.
    """
    found_segments = []
    for cyl in cyl_decomposition:
        regions_dict = cyl.regions
        veech_elem = cyl.veech_element
        for reg_1, reg_2 in itertools.product(regions_dict, repeat=2):
            # iterating over pairs of regions
            if debug:
                print(
                    f"Working on the possible region action"
                    f"{region_1} -> {region_2}")

            cylin_constraint = cyl.rat_ht_constraint
            # constraint coming from cylinder cyl
            other_constraint = regions_dict[reg_1].rat_ht_constraint
            # constraint from other cylinder that reg_1 is in
            acted_constraint = \
                regions_dict[reg_2].rat_ht_constraint.apply_matrix(veech_elem)
            # constraint from point after being acted on by veech_elem
            constraints = [cylin_constraint, other_constraint, acted_constraint]

#TODO: Add the constant offsets D:

            segment = produce_segment_three_polys(constraints)
            found_segments.append(segment)
            if debug:
                print(constraints, segment)









