from sage.all import *

from unittest import TestCase
from collections import defaultdict
from itertools import combinations

from samsurf import geom_equiv
from samsurf import comb_equiv
from samsurf import triangulation
from samsurf import unfolding


class Test(TestCase):
    def test_generate_code_marked(self):
        X = unfolding.triangulate_gothic1128(QQ(-1/2))
        codes_to_identifications = defaultdict(list)
        codes_to_vectors = defaultdict(list)
        for tri in range(len(X.triangles)):
            for edge in range(3):
                code, _ = geom_equiv.generate_code_marked(X, tri, edge)
                d = comb_equiv.canonical_relabel(X, tri, edge)
                dinv = {v: k for k, v in d.items()}

                ids = {d[e0]: d[e1] for e0, e1 in X.gluings.items()}
                codes_to_identifications[code].append(ids)

                vectors = []
                m = geom_equiv.get_normalization_matrix(X, tri, edge)
                for tri_idx in range(len(X.triangles)):
                    for edge_idx in range(3):
                        k0, k1 = dinv[(tri_idx, edge_idx)]
                        vectors.append(m*(X.triangles[k0][k1]))
                codes_to_vectors[code].append(vectors)
        for code, ids in codes_to_identifications.items():
            for ids0, ids1 in combinations(ids, 2):
                self.assertEqual(ids0, ids1)
        for code, vectors in codes_to_vectors.items():
            for vs0, vs1 in combinations(vectors, 2):
                self.assertEqual(vs0, vs1)
