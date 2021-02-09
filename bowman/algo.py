from collections import deque

from sage.all import *

from bowman.comb_equiv import gen_comb_equivs
from bowman.geom_equiv import gen_geom_equiv, gen_geom_equivs
from bowman.fund_dom import FundDom


def sigma(mat):
    [[a, b], [c, d]] = mat
    return sage.all.matrix([[a, -b], [-c, d]])


def _find_veech_equivs(idr0, idr1):
    if idr0.triangulation.code != idr1.triangulation.code:
        return []
    return gen_geom_equivs(idr0.triangulation, idr1.triangulation)


def generators_veech(trin):
    if not trin.is_delaunay():
        raise ValueError("Need to start with Delaunay Triangulation")
    r0 = trin.idr
    if r0.is_trivial:
        raise ValueError("The starting IDR is degenerate")

    if r0.has_self_equivalences:
        raise ValueError("Starting IDR has self-equivalences")

    edges_crossed = set()
    edges_zipped = {}
    idrs_to_visit = deque([r0])
    fund_dom = [r0]
    generators = []
    code_to_idr = {r0.triangulation.code: r0}

    while idrs_to_visit:
        idr_curr = idrs_to_visit.pop()
        for idx, edge in enumerate(idr_curr.polygon.edges):
            if edge not in edges_crossed and edge not in edges_zipped:
                edges_crossed |= {edge, edge.reverse()}
                idr_neighbor = idr_curr.cross_segment(idx)
                code_neighbor = idr_neighbor.triangulation.code
                if code_neighbor not in code_to_idr:
                    code_to_idr[code_neighbor] = idr_neighbor
                    fund_dom.append(idr_neighbor)
                    idrs_to_visit.appendleft(idr_neighbor)
                else:
                    t1 = idr_neighbor.triangulation
                    t2 = code_to_idr[code_neighbor].triangulation
                    ve = gen_geom_equiv(t1, t2)

                    m = sigma(ve)
                    edge_opp = edge.apply_mobius(m)
                    edges_zipped[edge] = edge_opp
                    edges_zipped[edge_opp] = edge
                    edges_zipped[edge.reverse()] = edge_opp.reverse()
                    edges_zipped[edge_opp.reverse()] = edge.reverse()
                    generators.append(m)

    return FundDom(fund_dom, generators, edges_zipped)
