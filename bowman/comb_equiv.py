import sage.all
from sage.all import *

import itertools
from collections import deque

from context import bowman
import bowman.triangulation
from bowman import triangulation


def homeo_respects_gluings(equiv, t1, t2):
    def image_edge(edge):
        perm = equiv[0]
        shift = equiv[1]
        return (perm[edge[0]], (edge[1] + shift[edge[0]]) % 3)

    # for all gluings (e1,e2) in t1, want to see if (F(e1), F(e2)) glued in t2
    for e1 in t1.edges:
        e2 = t1.gluings[e1]
        f1 = image_edge(e1)
        f2 = image_edge(e2)
        if t2.gluings[f1] != f2:
            return False

    return True


def gen_p_map_from_edge(t1, t2, f):
    def extend_map_to_nbrs(e, f):
        return {(e[0], (e[1] + k) % 3): (f[0], (f[1] + k) % 3) for k in range(3)}

    tris_to_visit = deque([0])
    partial_mapping = extend_map_to_nbrs((0, 0), f)
    # extend p_map to triangle 0, then mark it
    tris_visited = {0}

    while tris_to_visit:
        curr_tri = tris_to_visit.pop()
        if len(tris_visited) == len(t1.triangles):
            return partial_mapping
        for nbr_edge in [t1.gluings[(curr_tri, k)] for k in range(3)]:
            nbr_tri = nbr_edge[0]
            if nbr_tri not in tris_visited:
                im_nbr_edge = t2.gluings[(
                    partial_mapping[t1.gluings[(nbr_edge)]])]
                partial_mapping.update(
                    extend_map_to_nbrs(nbr_edge, im_nbr_edge))
                tris_visited.add(nbr_tri)
                tris_to_visit.appendleft(nbr_tri)

    return None


def gen_homeo_from_p_map(p_map, num_tris):
    # A P-map tells you, for each triangle, where *at least one* edge is mapped
    perm = []
    shifts = []

    for n in range(num_tris):
        edge = (0, 0)
        im_edge = (0, 0)
        # So: find an edge on tri n in the p_map
        for k in range(2):
            if (n, k) in p_map:
                edge = (n, k)
                im_edge = p_map[(n, k)]

        perm.append(im_edge[0])
        shifts.append((im_edge[1] - edge[1]) % 3)

    return (tuple(perm), tuple(shifts))


def gen_homeo_from_edge(t1, t2, e):
    partial_mapping = gen_p_map_from_edge(t1, t2, e)

    return gen_homeo_from_p_map(partial_mapping, len(t1.triangles))


# Returns a set of all combinatorial equivalences
def gen_comb_equivs(t1, t2):
    num_tris = len(t1.triangles)
    assert(num_tris == len(t2.triangles))

    poss_homeos = [gen_homeo_from_edge(t1, t2, e)
                   for e in itertools.product(range(num_tris), range(3))]

    return set(filter(lambda x: homeo_respects_gluings(x, t1, t2),
                      poss_homeos))

