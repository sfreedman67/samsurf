from sage.all import *

from unittest import TestCase
from collections import defaultdict
from itertools import combinations

from bowman import comb_equiv
from bowman import triangulation
from bowman import unfolding


class TestCombEquiv(TestCase):
    def test_canonical_relabel(self):
        octagon = triangulation.Triangulation.regular_octagon()
        answer = {(1, 1): (0, 0), (1, 2): (0, 1), (1, 0): (0, 2),
                  (0, 0): (1, 0), (0, 1): (1, 1), (0, 2): (1, 2),
                  (2, 2): (2, 0), (2, 0): (2, 1), (2, 1): (2, 2),
                  (4, 1): (3, 0), (4, 2): (3, 1), (4, 0): (3, 2),
                  (5, 0): (4, 0), (5, 1): (4, 1), (5, 2): (4, 2),
                  (3, 1): (5, 0), (3, 2): (5, 1), (3, 0): (5, 2)}
        self.assertEqual(comb_equiv.canonical_relabel(octagon, *(1, 1)), answer)

    def test_generate_code_marked(self):
        X = unfolding.triangulate_gothic1128(QQ(-2))
        codes_to_identifications = defaultdict(list)
        for tri in range(len(X.triangles)):
            for edge in range(3):
                code, _ = comb_equiv.generate_code_marked(X, tri, edge)
                d = comb_equiv.canonical_relabel(X, tri, edge)
                ids = {d[e0]: d[e1] for e0, e1 in X.gluings.items()}
                codes_to_identifications[code].append(ids)
        for code, ids in codes_to_identifications.items():
            for ids0, ids1 in combinations(ids, 2):
                self.assertEqual(ids0, ids1)

    def test_codes_comb(self):
        X = unfolding.triangulate_gothic1128(QQ(-1/2))
        codes_total = {comb_equiv.generate_code_marked(X, tri, edge)
                       for tri in range(len(X.triangles))
                       for edge in range(3)}
        assert X.codes_comb <= codes_total
