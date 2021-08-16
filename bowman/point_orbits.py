from sage.all import *

from bowman.triangle import Triangulation


def point_under_veech_action(trin, veech_elem, in_point):
    """
    Returns the location of a point under action of an element of veech group
    Inputs:
    * trin: Triangulation object
    * veech_elem: Matrix
    must be in the Veech group of trin
    * point: (int, 3-tuple)
    The first int represents the triangle id, the 3-tuple represents
    barycentric coordinates in that triangle.
    Output:
    out_point represented as (int, 3-tuple), where in_point ends up in trin
    under action of veech_elem
    """
    pass
    # Bug: Make Delaunay no longer preserves the marked points!! What do?


def points_preserved(trin, veech_elem, points_list)
# finds output of prev function for each point, iterating until the set is
# stable. Then we have all points preserved by a given veech_elem's subgroup

def reduce_periodic_points(trin, candidate_points)
# perform the above for every veech group generator
# output is list of periodic points!
