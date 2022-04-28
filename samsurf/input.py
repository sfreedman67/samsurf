import typing
import itertools

import sage.all
from sage.all import *
from scipy.spatial import Delaunay
import numpy as np

from samsurf import triangulation


def from_permutation(vectors: typing.Sequence, permutation: typing.Mapping[int, int]) -> triangulation.Triangulation:
    """Construct a translation surface via the permutation construction.

    Usage::

        >>> bottom_square = [sage.all.vector([1, 0]), sage.all.vector([0, 1])]
        >>> square = from_permutation(bottom_square, {0: 1, 1: 0})

    :param vectors: A list of 2-D Sage vectors forming one broken line.
    :param permutation: A dictionary from int to int representing a bijection on the set 1 to ``len(vectors)``.
    :rtype: A :class:`Triangulation` formed by the two broken lines
    """

    pts_bottom = list(itertools.accumulate(vectors))
    pts_top = list(itertools.accumulate(vectors[permutation[k]] for k in range(len(vectors))))

    # Don't double count the last vertex of each line
    pts = [sage.all.zero_vector(2), *pts_bottom, *reversed(pts_top[:-1])]

    trin = Delaunay(pts)

    tris = []
    for idx, (k0, k1, k2) in enumerate(trin.simplices):
        v0 = pts[k2] - pts[k1]
        v1 = pts[k0] - pts[k2]
        v2 = pts[k1] - pts[k0]
        tris.append(triangulation.Triangle(v0, v1, v2))

    gluings = {}
    vertex_to_pt_boundary = {}
    for s, t in itertools.product(range(len(tris)), range(3)):
        s1 = trin.neighbors[s][t]
        if s1 != -1:
            t1 = np.where(trin.neighbors[s1] == s)[0][0]
            gluings[(s, t)] = (s1, t1)
        else:
            vertex_to_pt_boundary[(s, t)] = tuple(pts[trin.simplices[s, (t + 1) % 3]])

    pt_boundary_to_vertex = {v: k for k, v in vertex_to_pt_boundary.items()}
    pt_opposite = {tuple(pt): tuple(pts[(k + len(pts) // 2) % len(pts)]) for k, pt in enumerate(pts)}
    for (s, t), pt in vertex_to_pt_boundary.items():
        gluings[(s, t)] = pt_boundary_to_vertex[pt_opposite[pt]]

    return triangulation.Triangulation(tris, gluings).make_nontrivial()


def regular_2kgon(k: int) -> triangulation.Triangulation:
    """Triangulate a regular 2k-gon"""
    zetaN = QQbar(CyclotomicField(2 * k).gen())
    vectors = []
    for x in range(k):
        w = zetaN ** x
        u = (w + w.conjugate()) / QQbar(2)
        v = (w - w.conjugate()) / QQbar(2 * i)
        vectors.append(sage.all.vector([AA(u), AA(v)]))
    return from_permutation(vectors, {x: (k - 1 - x) for x in range(k)})


if __name__ == "__main__":
    # TODO: Cleaning, big time
    X = regular_2kgon(k=5)

    # might have to try multiple shears
    shear = sage.all.matrix([[1, QQ(1/5)], [0, 1]])
    X = X.apply_matrix(shear).make_delaunay().apply_matrix(shear.inverse())

    from comb_equiv import canonical_relabel
    cr = canonical_relabel(X, 0, 0)
    perm, shifts = zip(*[cr[(s, 0)] for s in range(len(X.triangles))])

    tris_shifted = []
    for tri, shift in zip(X.triangles, shifts):
        vectors_shifted = tri[-shift:] + tri[:-shift]
        tri_shifted = triangulation.Triangle(*vectors_shifted)
        tris_shifted.append(tri_shifted)

    tris_relabelled = [tris_shifted[perm.index(k)] for k in range(len(perm))]

    gluings_relabelled = {cr[x]: cr[y] for x, y in X.gluings.items()}
    Y = triangulation.Triangulation(tris_relabelled, gluings_relabelled)
