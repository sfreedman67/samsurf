import sage.all
from sage.all import *

import itertools
from unittest import TestCase

from bowman import comb_equiv
from bowman import triangulation


class TestCombEquiv(TestCase):

    # TODO take triangulated surface, then relabel everything

    def test_gen_homeo_from_edge(self):
        sq_torus = triangulation.Triangulation.square_torus()

        comb_equiv1 = comb_equiv.CombEquiv.from_edge(sq_torus, sq_torus, (0, 0))
        matching2 = comb_equiv.CombEquiv.from_edge(sq_torus, sq_torus, (1, 1))

        self.assertEqual(comb_equiv1[:2], ((0, 1), (0, 0)))
        self.assertEqual(matching2[:2], ((1, 0), (1, 1)))

    def test_torus_self_equivs(self):
        sq_torus = triangulation.Triangulation.square_torus()
        answer = list(itertools.product([(1, 0), (0, 1)],
                                        [(i, i) for i in range(3)]))
        output = [x[:2] for x in
                  comb_equiv.gen_comb_equivs(sq_torus, sq_torus)]
        self.assertCountEqual(output, answer)

    def test_veech_octagon_no_comb_autos(self):
        octagon = triangulation.Triangulation.regular_octagon()
        answer = comb_equiv.gen_comb_equivs(octagon, octagon)[0][:2]
        self.assertEqual(answer,
                         ((0, 1, 2, 3, 4, 5), (0, 0, 0, 0, 0, 0)))

    def test_deduce_rotation(self):
        tri_equil = triangulation.Triangle(vector([1, 0]), vector([-1, 1]), vector([0, -1]))
        tris = [tri_equil] * 4
        gluings = {(0, 0): (3, 0), (0, 1): (1, 0), (0, 2): (2, 0),
                   (1, 0): (0, 1), (1, 1): (2, 2), (1, 2): (3, 1),
                   (2, 0): (0, 2), (2, 1): (3, 2), (2, 2): (1, 1),
                   (3, 0): (0, 0), (3, 1): (1, 2), (3, 2): (2, 1)}
        t = triangulation.Triangulation(tris, gluings, QQ)

        assert False, "Not a test"

