import sage.all
from sage.all import *


def gen_geom_equiv(trin1, trin2, ce):
    m = _get_matrix_from_tris(trin1.triangles[0],
                              trin2.triangles[ce.perm[0]],
                              ce.shift[0])

    return m if _agrees_on_tris(m, ce) else None


def _get_matrix_from_tris(tri1, tri2, shift):
    vs = (tri1[0], -tri1[2])
    ws = (tri2[shift], -tri2[(shift + 2) % 3])
    return _get_matrix_from_vectors(vs, ws)


def _get_matrix_from_vectors(pair1, pair2):
    """get matrix sending (v1, v2) to (w1, w2)"""
    (a1, a2), (b1, b2) = pair1
    (c1, c2), (d1, d2) = pair2

    m1 = 1/(a1 * b2 - a2 * b1) * sage.all.matrix([[b2, -b1], [-a2, a1]])
    m2 = sage.all.matrix([[c1, d1], [c2, d2]])

    return m2 * m1


def _agrees_on_tris(m, ce):
    return all(_agrees_on_tri(m, ce, idx)
               for idx in range(len(ce.source.triangles)))


def _agrees_on_tri(m, ce, idx):
    tri1 = ce.source.triangles[idx]
    perm, shift = ce.perm[idx], ce.shift[idx]
    tri2 = ce.target.triangles[perm]

    v1, v2 = tri1[0], -tri1[2]
    w1, w2 = tri2[shift], -tri2[(shift + 2) % 3]
    return m * v1 == w1 and m * v2 == w2
