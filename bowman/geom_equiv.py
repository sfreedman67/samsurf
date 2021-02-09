import sage.all
from sage.all import *

from bowman.comb_equiv import canonical_relabel


def gen_geom_equiv(trin1, trin2):
    if trin1.code_geom != trin2.code_geom:
        return None

    _, (tri1, edge1) = next(iter(trin1.codes_geom))
    _, (tri2, edge2) = next(iter(trin2.codes_geom))

    m1 = get_normalization_matrix(trin1, tri1, edge1)
    m2 = get_normalization_matrix(trin2, tri2, edge2)

    return m2.inverse() * m1


def gen_geom_equivs(trin1, trin2):
    if trin1.code_geom != trin2.code_geom:
        return []

    _, (tri1, edge1) = next(iter(trin1.codes_geom))
    m1 = get_normalization_matrix(trin1, tri1, edge1)

    ges = []
    for _, (tri2, edge2) in trin2.codes_geom:
        m2 = get_normalization_matrix(trin2, tri2, edge2)
        ge = m2.inverse() * m1
        if -ge not in ges:
            ges.append(ge)

    return ges


def get_normalization_matrix(t, tri, edge):
    v1 = t.triangles[tri][edge]
    v2 = -t.triangles[tri][(edge + 2) % 3]
    return sage.all.matrix([v1, v2]).transpose().inverse()


def generate_code_marked(t, tri, edge):
    relabel = canonical_relabel(t, tri, edge)
    relabel_inv = {v: k for k, v in relabel.items()}

    m = get_normalization_matrix(t, tri, edge)
    frames = []
    for k in range(len(t.triangles)):
        tri1, edge1 = relabel_inv[(k, 0)]
        v = t.triangles[tri1][edge1]
        w = -t.triangles[tri1][(edge1 + 2) % 3]
        frames.append(sage.all.matrix([v, w]).transpose())

    frames_new = [m * frame for frame in frames]
    for frame in frames_new:
        frame.set_immutable()

    code = tuple(frame for frame in frames_new)

    return hash(code), (tri, edge)
