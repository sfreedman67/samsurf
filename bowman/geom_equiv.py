import sage.all
from sage.all import *

import bowman.comb_equiv as ce


def gen_geom_equivs(trin1, trin2):
    """
    A geometrix equivalence is a matrix as well as a dictionary: edge -> edge
    The matrix is applied to all triangles in trin1, to produce a triangulation
    which is exactly trin2, upto relabeling.
    The output dictionary tells us how to relabel
    """
    if trin1.code != trin2.code:
        return []

    _, (tri1, edge1) = next(iter(trin1.codes_geom))
    m1 = get_normalization_matrix(trin1, tri1, edge1)
    # matrix sending trin1 to its min geometric equivalence
    rel1 = ce.canonical_relabel(trin1, tri1, edge1)
    # relabeling sending the edges of trin1 to the min geometric equivalence

    ges = []
    for _, (tri2, edge2) in trin2.codes_geom:
        m2 = get_normalization_matrix(trin2, tri2, edge2)
        # matrix sending trin2 to its min geometric equivalence
        matrix_composed = m2.inverse() * m1
        # matrix sending trin1 to trin2

        rel2 = ce.canonical_relabel(trin2, tri2, edge2)
        # relabeling sending the edges of trin2 to min geometric equivalence
        rel2_inv = {val: key for key, val in rel2.items()}
        relabel_composed = {key: rel2_inv[rel1[key]] for key in rel1}
        # relabeling sending trin1 edges to trin2 edges

        ge = matrix_composed, relabel_composed  # geometric equivalence
        ges.append(ge)
    return ges


def gen_geom_equiv(trin1, trin2):
    equivs = gen_geom_equivs(trin1, trin2)
    return None if not equivs else equivs[0]


def get_normalization_matrix(t, tri, edge):
    v1 = t.triangles[tri][edge]
    v2 = -t.triangles[tri][(edge + 2) % 3]
    return sage.all.matrix([v1, v2]).transpose().inverse()


def generate_code_marked(t, tri, edge):
    relabel = ce.canonical_relabel(t, tri, edge)
    relabel_inv = {v: k for k, v in relabel.items()}

    frames = []
    for k in range(len(t.triangles)):
        tri1, edge1 = relabel_inv[(k, 0)]
        v = t.triangles[tri1][edge1]
        w = -t.triangles[tri1][(edge1 + 2) % 3]
        frames.append(sage.all.matrix([v, w]).transpose())

    m = get_normalization_matrix(t, tri, edge)
    frames_new = [m * frame for frame in frames]
    for frame in frames_new:
        frame.set_immutable()

    code = tuple(frame for frame in frames_new)
    return hash(code), (tri, edge)
