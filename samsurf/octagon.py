import collections

import sage.all
from sage.all import *

from bowman import triangulation

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


def octagon(a, b):
    e1, e2 = sage.all.vector([1, 0]), sage.all.vector([0, 1])

    triangles = [triangulation.Triangle(e1, e2, -e1 - e2),
                 triangulation.Triangle(-e1, -e2, e1 + e2),
                 triangulation.Triangle(a, -a - e1 - e2, e1 + e2),
                 triangulation.Triangle(-a, a + e1 + e2, -e1 - e2),
                 triangulation.Triangle(-b, e1 + e2 + a, -a - e2 - e1 + b),
                 triangulation.Triangle(b, -e1 - e2 - a, -b + e1 + e2 + a)]

    gluings = {(0, 0): (1, 0), (0, 1): (1, 1), (0, 2): (2, 2), (1, 2): (3, 2),
               (2, 0): (3, 0), (2, 1): (4, 1), (3, 1): (5, 1), (4, 0): (5, 0), (4, 2): (5, 2)}
    gluings.update({v: k for k, v in gluings.items()})

    return triangulation.Triangulation(triangles, gluings)


if __name__ == "__main__":
    ax, ay = var("ax, ay")
    bx, by = var("bx, by")
    a = sage.all.vector([ax, ay])
    b = sage.all.vector([bx, by])
    X = octagon(a, b)

    equations = [tuple(x.full_simplify() if x not in QQ else x for x in hinge._coefficients)
                 for hinge in X.hinges]

    eqs = [((a + c).full_simplify() if a + c not in QQ else a + c) == 0 for a, b, c in equations]
    print(solve(eqs[2:], [ax, ay, bx, by]))





