import sage.all
from sage.all import *

import unittest
import itertools

from context import bowman
import bowman.triangulation
from bowman import triangulation
import bowman.comb_equiv
from bowman import comb_equiv

class CombEquivTests(unittest.TestCase):

    # TODO take triangulated surface, then relabel everything

    def test_gen_homeo_from_edge(self):
        sq_torus = triangulation.Triangulation.square_torus()

        homeo1 = comb_equiv.gen_homeo_from_edge(sq_torus, sq_torus, (0, 0))
        homeo2 = comb_equiv.gen_homeo_from_edge(sq_torus, sq_torus, (1, 1))

        self.assertEqual(homeo1, ((0, 1), (0, 0)))
        self.assertEqual(homeo2, ((1, 0), (1, 1)))

    def test_torus_self_equivs(self):
        sq_torus = triangulation.Triangulation.square_torus()
        equivs_two_tris = set(itertools.product([(1, 0), (0, 1)], [(i, i)
                                                                   for i in range(3)]))
        self.assertEqual(comb_equiv.gen_comb_equivs(
            sq_torus, sq_torus), equivs_two_tris)

    def test_veech_octagon_no_comb_autos(self):
        octagon = triangulation.Triangulation.regular_octagon()
        comb_equiv.gen_comb_equivs(octagon, octagon)
        self.assertEqual(comb_equiv.gen_comb_equivs(octagon, octagon),
                         {((0, 1, 2, 3, 4, 5), (0, 0, 0, 0, 0, 0))})


if __name__ == "__main__":
    unittest.main(verbosity=2)

    octagon_suite = unittest.TestSuite()
    octagon_suite.addTest(MyTests("test_veech_octagon_no_comb_autos"))
    runner = unittest.TextTestRunner()
    # runner.run(octagon_suite)