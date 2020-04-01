# N = number of triangles , numbered 0 \le i < N
# j = 0,1,2 denotes side of triangle
# Must have 2E = 3T

# A combinatorial equivalence consists of a bijection on the triangles that respects edge gluings
# When we biject triangles, we need a "shift" to tell how we glued


import itertools as it
import unittest
from sage.all import *
import flatsurf
import sage.graphs.generic_graph as sg
from collections import deque


def homeo_respects_gluings(equiv, t1, t2):
    # equiv[0] = perm, equiv[1] = shift
    perm = equiv[0]
    shift = equiv[1]

    def image_edge(edge):
        return (perm[edge[0]], (edge[1] + shift[edge[0]]) % 3)

    # for all gluings (e1,e2) in t1, want to see if (F(e1), F(e2)) glued in t2
    for e1 in t1.edge_iterator():
        e2 = t1.opposite_edge(e1)
        f1 = image_edge(e1)
        f2 = image_edge(e2)
        if t2.opposite_edge(f1) != f2:
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
        if len(tris_visited) == t1.num_polygons():
            return partial_mapping
        for nbr_edge in [t1.opposite_edge(curr_tri, k) for k in range(3)]:
            nbr_tri = nbr_edge[0]
            if nbr_tri not in tris_visited:
                im_nbr_edge = t2.opposite_edge(
                    partial_mapping[t1.opposite_edge(nbr_edge)])
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

    return gen_homeo_from_p_map(partial_mapping, t1.num_polygons())


# Returns a list of all combinatorial equivalences
def gen_comb_equivs(t1, t2):
    num_tris = t1.num_polygons()
    assert(num_tris == t2.num_polygons())

    poss_homeos = [gen_homeo_from_edge(t1, t2, e)
                   for e in it.product(range(num_tris), range(3))]

    return set(filter(lambda x: homeo_respects_gluings(x, t1, t2),
                      poss_homeos))


class MyTests(unittest.TestCase):

    def test_gen_homeo_from_edge(self):
        sq_torus = flatsurf.translation_surfaces.square_torus()
        t = sq_torus.delaunay_triangulation()

        homeo1 = gen_homeo_from_edge(t, t, (0, 0))
        homeo2 = gen_homeo_from_edge(t, t, (1, 1))

        self.assertEqual(homeo1, ((0, 1), (0, 0)))
        self.assertEqual(homeo2, ((1, 0), (1, 1)))

    def test_torus_self_equivs(self):
        t = flatsurf.translation_surfaces.square_torus()
        DT_sq_torus = t.delaunay_triangulation()
        equivs_two_tris = set(it.product([(1, 0), (0, 1)], [(i, i)
                                                            for i in range(3)]))
        self.assertEqual(gen_comb_equivs(
            DT_sq_torus, DT_sq_torus), equivs_two_tris)

    def test_veech_octagon_no_comb_autos(self):
        v_oct = flatsurf.translation_surfaces.regular_octagon()
        DT_v_oct = v_oct.delaunay_triangulation()
        gen_comb_equivs(DT_v_oct, DT_v_oct)
        self.assertEqual(gen_comb_equivs(DT_v_oct, DT_v_oct),
                         {((0, 1, 2, 3, 4, 5), (0, 0, 0, 0, 0, 0))})


unittest.main(verbosity=2)

octagon_suite = unittest.TestSuite()
octagon_suite.addTest(MyTests("test_veech_octagon_no_comb_autos"))
runner = unittest.TextTestRunner()
# runner.run(octagon_suite)
