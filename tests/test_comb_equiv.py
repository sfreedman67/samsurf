import sage.all
from sage.all import *

import itertools
from unittest import TestCase

from bowman import comb_equiv
from bowman.comb_equiv import *
from bowman import triangulation


class TestCombEquiv(TestCase):

    def test_gen_homeo_from_edge(self):
        sq_torus = triangulation.Triangulation.square_torus()

        ce1 = comb_equiv.CombEquiv.from_edge(sq_torus, sq_torus, (0, 0))
        ce2 = comb_equiv.CombEquiv.from_edge(sq_torus, sq_torus, (1, 1))

        self.assertEqual(ce1[:2], ((0, 1), (0, 0)))
        self.assertEqual(ce2[:2], ((1, 0), (1, 1)))

    def test_torus_self_equivs(self):
        sq_torus = triangulation.Triangulation.square_torus()
        answer = list(itertools.product([(1, 0), (0, 1)],
                                        [(i, i) for i in range(3)]))
        output = [x[:2] for x in
                  comb_equiv.gen_comb_equivs(sq_torus, sq_torus)]
        self.assertCountEqual(output, answer)

    def test_octagon_self_equivs(self):
        octagon = triangulation.Triangulation.regular_octagon()
        identity, rotate = comb_equiv.gen_comb_equivs(octagon, octagon)
        self.assertEqual(identity[:2], ((0, 1, 2, 3, 4, 5), (0, 0, 0, 0, 0, 0)))
        self.assertEqual(rotate[:2], ((1, 0, 3, 2, 5, 4), (0, 0, 0, 0, 0, 0)))

    def test_canonical_relabel(self):
        octagon = triangulation.Triangulation.regular_octagon()
        answer = {(1, 1): (0, 0), (1, 2): (0, 1), (1, 0): (0, 2),
                  (0, 1): (1, 0), (0, 2): (1, 1), (0, 0): (1, 2),
                  (3, 2): (2, 1), (3, 0): (2, 2), (3, 1): (2, 0),
                  (2, 2): (3, 1), (2, 0): (3, 2), (2, 1): (3, 0),
                  (5, 1): (4, 0), (5, 2): (4, 1), (5, 0): (4, 2),
                  (4, 1): (5, 0), (4, 2): (5, 1), (4, 0): (5, 2)}
        self.assertEqual(comb_equiv.canonical_relabel(octagon, *(1, 1)), answer)

    def test_generate_code(self):
        octagon = triangulation.Triangulation.regular_octagon()
        octagon_flipped = octagon.flip_hinge((0, 2))

        self.assertFalse(gen_comb_equivs(octagon, octagon_flipped))
        for tri, edge in itertools.product(range(6), range(3)):
            self.assertNotEqual(generate_code_marked(octagon, 0, 0),
                                generate_code_marked(octagon_flipped, tri, edge))

        self.assertNotEqual(octagon.code_comb, octagon_flipped.code_comb)

        mc_l = triangulation.Triangulation.mcmullen_l(QQ(3), QQ(3))
        idr = mc_l.idr
        t = idr.triangulation
        idr0 = idr.cross_segment(0)
        t0 = idr0.triangulation

        for tri, edge in itertools.product(range(10), range(3)):
            if (tri, edge) != (7, 0) and (tri, edge) != (2, 0):
                self.assertNotEqual(generate_code_marked(t, 5, 0),
                                    generate_code_marked(t0, tri, edge))
            else:
                self.assertEqual(generate_code_marked(t, 5, 0),
                                 generate_code_marked(t0, tri, edge))

        self.assertEqual(t.code_comb, t0.code_comb)
