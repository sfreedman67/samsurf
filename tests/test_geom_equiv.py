from unittest import TestCase

import sage.all
from sage.all import *

from context import bowman
import bowman.geom_equiv
from bowman import geom_equiv


class Test(TestCase):
    def test_get_matrix(self):
        v1 = vector([2, 1])
        v2 = vector([1, 1])

        w1 = v2
        w2 = vector([1, 2])

        ans = matrix([[0, 1], [-1, 3]])
        self.assertEqual(geom_equiv._get_matrix_from_vectors((v1, v2), (w1, w2)), ans)
