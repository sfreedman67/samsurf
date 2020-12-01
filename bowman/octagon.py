import collections

import sage.all
from sage.all import *
from sage.groups.matrix_gps.finitely_generated import MatrixGroup

from bowman import triangulation
from bowman import algo

K = QuadraticField(2)
A = K.gen()

VERTS = [
    sage.all.vector([1, 0]),
    sage.all.vector([A / 2, A / 2]),
    sage.all.vector([0, 1]),
    sage.all.vector([-A / 2, A / 2]),
    sage.all.vector([-1, 0]),
    sage.all.vector([-A / 2, -A / 2]),
    sage.all.vector([0, -1]),
    sage.all.vector([A / 2, -A / 2])
]

# noinspection PyArgumentList
PC = PointConfiguration(VERTS, connected=True, fine=True, star=None, regular=None)

Trin = collections.namedtuple("Trin", ["comb", "geom"])


def tri_geom(tri_comb):
    n0, n1, n2 = tri_comb
    return triangulation.Triangle(VERTS[n1] - VERTS[n0],
                                  VERTS[n2] - VERTS[n1],
                                  VERTS[n0] - VERTS[n2])


def opposite(edge_comb):
    start, stop = edge_comb
    if (stop - start) % 8 == 1:
        return (start + 4) % 8, (stop + 4) % 8
    else:
        return stop, start


def edges_comb(tri_comb):
    n0, n1, n2 = tri_comb
    return [(n0, n1), (n1, n2), (n2, n0)]


def gluings(trin_comb):
    comb_to_lab = {edge_comb: (tri_lab, edge_lab)
                   for tri_lab, tri_comb in enumerate(trin_comb)
                   for edge_lab, edge_comb in enumerate(edges_comb(tri_comb))}

    return {(tri_lab, edge_lab): comb_to_lab[opposite(edge_comb)]
            for tri_lab, tri_comb in enumerate(trin_comb)
            for edge_lab, edge_comb in enumerate(edges_comb(tri_comb))}


def geom_from_comb(trin_comb):
    tris_geom = [tri_geom(tri_comb) for tri_comb in trin_comb]

    return triangulation.Triangulation(tris_geom, gluings(trin_comb))


def gen_trins():
    return [Trin(list(trin_comb), geom_from_comb(list(trin_comb)))
            for trin_comb in PC.triangulations_list()]


def get_idr(trin):
    return trin.geom.idr


def gen_idrs_nontrivial():
    idrs = [get_idr(trin) for trin in gen_trins()]
    return [IDR for IDR in idrs if not IDR.is_trivial]


def gen_trins_with_idr_nontrivial():
    return [trin for trin in gen_trins() if not get_idr(trin).is_trivial]


def plot_trin_comb(trin_comb):
    octagon = sage.all.plot.polygon.polygon(VERTS, fill=False)

    for tri_comb in trin_comb:
        for start, stop in edges_comb(tri_comb):
            if (stop - start) % 8 != 1:
                octagon += sage.all.plot.line.line([VERTS[start], VERTS[stop]])

    return octagon


if __name__ == "__main__":
    X = gen_trins_with_idr_nontrivial()[0].geom
    r0 = X.idr
    r1, r2, r3 = r0.neighbors

    # x0, y0 = var('x0, y0')
    # x1, y1 = var('x1, y1')
    # x2, y2 = var('x2, y2')
    # x3, y3 = var('x3, y3')
    #
    # z0 = vector([x0, y0])
    # z1 = vector([x1, y1])
    # z2 = vector([x2, y2])
    # z3 = vector([x3, y3])
    #
    # triangles = [triangulation.Triangle(z0, z1, -z1 - z0),
    #              triangulation.Triangle(z2, -z0 - z1 - z2, z0 + z1),
    #              triangulation.Triangle(-z3 - z2 - z1 - z0, z0 + z1 + z2, z3),
    #              triangulation.Triangle(z0 + z1 + z2 + z3, -z0, -z3 - z2 - z1),
    #              triangulation.Triangle(-z1, -z2 - z3, z1 + z2 + z3),
    #              triangulation.Triangle(-z3, z2 + z3, -z2)]
    #
    # gluings = {(0, 0): (3, 1), (0, 1): (4, 0), (0, 2): (1, 2),
    #            (1, 0): (5, 2), (1, 1): (2, 1), (2, 0): (3, 0),
    #            (2, 2): (5, 0), (3, 2): (4, 2), (4, 1): (5, 1)}
    # gluings.update({v: k for k, v in gluings.items()})
    #
    # X = triangulation.Triangulation(triangles, gluings, sage.all.SR)
    # hinge0 = X.hinges[0]
    # a, b, c = hinge0.halfplane
    # print(hinge0.incircle_det.full_simplify())


