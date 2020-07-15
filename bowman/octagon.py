import sage.all
from sage.all import *

import functools
import collections

import triangulation
from triangulation import Triangle, Triangulation

import intersect_halfplanes
from intersect_halfplanes import intersect_halfplanes

from polygon import plot_polygon

K = QuadraticField(2)
A = K.gen()

VERTS = [
    vector([1, 0]),
    vector([A / 2, A / 2]),
    vector([0, 1]),
    vector([-A / 2, A / 2]),
    vector([-1, 0]),
    vector([-A / 2, -A / 2]),
    vector([0, -1]),
    vector([A / 2, -A / 2])
]

PC = PointConfiguration(VERTS)

Trin = collections.namedtuple("Trin", ["comb", "geom"])


def tri_geom(tri_comb):
    n0, n1, n2 = tri_comb
    return Triangle((VERTS[n1] - VERTS[n0],
                     VERTS[n2] - VERTS[n1],
                     VERTS[n0] - VERTS[n2]))


def opposite(edge_comb):
    start, end = edge_comb
    if (end - start) % 8 == 1:
        return ((start + 4) % 8, (end + 4) % 8)
    else:
        return (end, start)


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

    return Triangulation(tris_geom, gluings(trin_comb), K)


def gen_trins():
    return [Trin(list(trin_comb), geom_from_comb(list(trin_comb)))
            for trin_comb in PC.triangulations_list()]


def IDR(trin):
    return intersect_halfplanes(trin.geom.halfplanes())


def is_nontriv(IDR):
    return IDR is not None and len(IDR) > 2


def gen_IDRs_nontriv():
    IDRs = [IDR(trin) for trin in gen_trins()]
    return [IDR for IDR in IDRs if is_nontriv(IDR)]


def gen_trins_IDR_nontriv():
    return [trin for trin in gen_trins() if is_nontriv(IDR(trin))]


def plot_trin_comb(trin_comb):
    octagon = sage.plot.polygon.polygon(VERTS, fill=False)

    for tri_comb in trin_comb:
        for start, end in edges_comb(tri_comb):
            if (end - start) % 8 != 1:
                octagon += sage.plot.line.line([VERTS[start], VERTS[end]])

    return octagon


if __name__ == "__main__":
    _, trin_geom = gen_trins_IDR_nontriv()[0]
    for halfplane, hinges in trin_geom.halfplanes_to_hinges.items():
        print(halfplane, hinges)

