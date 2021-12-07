# Place to store the prototype of the Rational Height Lemma

import itertools
from sage.all import *
from samsurf.linear_xy_poly import LinearXYPoly
from samsurf.cylinder import Cylinder
from sage.rings.number_field.number_field import is_QuadraticField


def common_field(list_of_fields):
    # finds the largest field in list_of_fields
    #TODO: only works with quadratic fields, for now
    output = QQ
    found_quadratic = False
    for field in list_of_fields:
        if field == QQ:
            continue
        else:
            if is_QuadraticField(field):
                if not found_quadratic:
                    output = field  # update when first quadratic field found
                    found_quadratic = True
                elif field.discriminant() == output.discriminant():
                    continue
                else:
                    raise ValueError(f"Input {field} is different field from {output}")
            else:
                raise ValueError(f"Inputs must be quadratic field, which {field} is not.")
    return output


def euclidean_gcd(a, b):
    # returns the gcd of number field elements a, b using the Euclidean algorithm
    if b > a:
        a, b = b, a
    if b == 0:
        return a
    else:
        return euclidean_gcd(b, a - b)


def list_lcm(list_of_nums):
    """
    Returns the lcm of a list of numbers
    """
    list_so_far = list_of_nums[:]  # deep copy of the list
    lcm = list_so_far.pop()
    while list_so_far:
        new_num = list_so_far.pop()
        lcm = lcm * new_num / euclidean_gcd(lcm, new_num)
    return lcm



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

    base_field = common_field(poly.base_field for poly in list_of_poly)
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


def reduce_cylinder_constraints(cylinder, global_veech_elem, debug=False):
    """
    Finds all the segments on a cylinder where periodic points might lie.
    Inputs:
    * cylinder: Cylinder object
        Represents all the cylinders in a given direction
    global_veech_elem is the parabolic element along the cylinder direction
    for the whole surface.
    Output:
        List of line segments where periodic points lie.
        Two segments are produced for each pair of regions in a given cylinder
    #TODO: This code assumes cylinder direction is horizontal, and other
    direction is vertical.
    """

    found_segments = {}
    farthest_idx = ((cylinder.num_twists(global_veech_elem) + 1) *
                    cylinder.num_triangles)
    # farthest triangle a point travels to under action of global_veech_elem
    if debug:
        print(f"----------------")
        print(f"Looking at cylinder {list(cylinder.idx_dict.keys())} with")
        print(f"direction {cylinder.direction}")
        print(f"and veech elem {global_veech_elem}\n")
        print(f"Points in cylinder travel upto {farthest_idx} triangles away\n")

    for idx, tri in cylinder.idx_dict.items():
        tri_segments = []
        cylin_constraint = tri.constraint_in_direction(cylinder.direction)
        other_constraint = tri.constraint_in_direction(cylinder.other_direction)
        if debug:
            print(f"--------")
            print(f"working on triangle {idx}")
            print(f"got constraints {cylin_constraint} and {other_constraint}")

        for i in range(0, int(farthest_idx), 2):
            final_constraint = cylinder.get_third_constraint(tri, i, global_veech_elem)
            segment = produce_segment_three_polys(
                [cylin_constraint, other_constraint, final_constraint])
            if debug:
                print(f"Considering the triangle {i} indices away.")
                print(f"Final constraint is {final_constraint}.")
                print(f"The case {idx} -> {idx + i} produces line {segment}\n")
            tri_segments.append(segment)

        found_segments[idx] = tri_segments

    return found_segments


def bicuspid_segments(triangulation, debug=False):
    """
    Given a bicuspid triangulation with both horizontal and vertical cylinders
    finds the segments of candidate periodic points in each triangle.
    """

    hor_dir = vector([1, 0])
    ver_dir = vector([0, 1])

    refined_tris, ver_cyl_idxs = \
        triangulation.make_directional_triangulation(ver_dir)
    assert refined_tris.check_horiz(), "Triangulation not bicuspid"
    #TODO Assumes that the vertical sheared triangulation will have horizontal
    # cylinders too - might be wrong
    hor_cyl_idxs = refined_tris.get_horiz_cyls()

    hor_cyls = [Cylinder.from_indices(hor_dir, refined_tris, cyl)
                for cyl in hor_cyl_idxs]
    hor_modulus_lcm = list_lcm([1 / cyl.modulus for cyl in hor_cyls])
    hor_veech_elem = matrix([[1, hor_modulus_lcm], [0, 1]])
    hor_segments = {}
    for cyl in hor_cyls:
        hor_segments.update(
            reduce_cylinder_constraints(cyl, hor_veech_elem, debug=debug))

    ver_cyls = [Cylinder.from_indices(ver_dir, refined_tris, cyl)
                for cyl in ver_cyl_idxs]
    ver_modulus_lcm = list_lcm([1 / cyl.modulus for cyl in ver_cyls])
    ver_veech_elem = matrix([[1, 0], [ver_modulus_lcm, 1]])
    ver_segments = {}
    for cyl in ver_cyls:
        ver_segments.update(
            reduce_cylinder_constraints(cyl, ver_veech_elem, debug=debug))

    # if debug:
    #     print(f"Horizontal veech element used")

    tri_segments = {idx: (hor_segments[idx], ver_segments[idx])
                    for idx in hor_segments}
    return tri_segments


def segments_for_plotting(tri_segments):
    """
    Takes the output of bicuspid_segments, and produces a simple
    dictionary which can be plotted.
    Loses information for triangles where both directions yield valid lists
    #TODO: once intersection of marked lines implemented, don't lose the info
    """
    output = {}
    for idx, (hor_segments, ver_segments) in tri_segments.items():
        if 0 in hor_segments:
            if 0 in ver_segments:
                raise ValueError("Triangle {idx} has zero in both lists")
            else:
                tri_list = ver_segments[:]
        else:
            if 0 in ver_segments:
                tri_list = hor_segments[:]
            else:
                tri_list = (ver_segments[:]
                            if len(hor_segments) > len(ver_segments)
                            else hor_segments[:])
        output[idx] = [LinearXYPoly.from_polynomial(poly).get_coeffs()
                       for poly in tri_list]
    return output
