from unittest import TestCase

from sage.all import *

from bowman.mobius import apply_mobius
from bowman.polygon import Point
from bowman.radical import Radical


class Test(TestCase):
    def test_apply_mobius(self):
        m1 = sage.all.matrix([[1, 7], [0, 1]])
        m2 = sage.all.matrix([[2, 1], [1, 1]])

        self.assertEqual(apply_mobius(m1, oo), oo)
        self.assertEqual(apply_mobius(m1, Point(3, 4)), Point(10, 4))
        self.assertEqual(apply_mobius(m1, Point(Radical(2, 1, 3), 0)),
                         Point(Radical(9, 1, 3), 0))

        self.assertEqual(apply_mobius(m2, oo), Point(2, 0))
        self.assertEqual(apply_mobius(m2, Point(3, 4)),
                         Point(QQ(36 / 20), QQ(1 / 100)))
        self.assertEqual(apply_mobius(m2, Point(-1, 0)), oo)
        self.assertEqual(apply_mobius(m2, Point(Radical(2, -1, 3), 0)),
                         Point(Radical(QQ(3 / 2), QQ(-1 / 6), 3), 0))
