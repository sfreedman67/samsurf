from collections import deque

from sage.all import *

from bowman.comb_equiv import gen_comb_equivs, canonical_relabel
from bowman.geom_equiv import gen_geom_equiv
from bowman.fund_dom import FundDom


def sigma(mat):
    [[a, b], [c, d]] = mat
    return sage.all.matrix([[a, -b], [-c, d]])


def _find_veech_equiv(idr0, idrs):
    for idr in idrs:
        t0, t = idr0.triangulation, idr.triangulation
        if t0.code_comb == t.code_comb and t0.code_geom == t.code_geom:
            return gen_geom_equiv(t0, t)
    return None


def _find_veech_equivs(idr0, idr1):
    ces = gen_comb_equivs(idr0.triangulation, idr1.triangulation)
    ges = []
    if ces:
        ge = gen_geom_equiv(idr0.triangulation, idr1.triangulation)
        if ge is not None and -ge not in ges:
            ges.append(ge)
    return ges


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
    gens = []
    codes_to_idrs = {r0.triangulation.code_comb: [r0]}

    while idrs_to_visit:
        idr_curr = idrs_to_visit.pop()
        for idx, edge in enumerate(idr_curr.polygon.edges):
            if edge not in edges_crossed and edge not in edges_zipped:
                edges_crossed |= {edge, edge.reverse()}
                idr_neighbor = idr_curr.cross_segment(idx)

                code_neighbor = idr_neighbor.triangulation.code_comb
                if code_neighbor not in codes_to_idrs:
                    codes_to_idrs[code_neighbor] = [idr_neighbor]
                    fund_dom.append(idr_neighbor)
                    idrs_to_visit.appendleft(idr_neighbor)
                else:
                    ve = _find_veech_equiv(idr_neighbor, codes_to_idrs[code_neighbor])
                    if ve is not None:
                        m = sigma(ve)
                        edge_opp = edge.apply_mobius(m)
                        edges_zipped[edge] = edge_opp
                        edges_zipped[edge_opp] = edge
                        edges_zipped[edge.reverse()] = edge_opp.reverse()
                        edges_zipped[edge_opp.reverse()] = edge.reverse()
                        gens.append(m)
                    else:
                        codes_to_idrs[code_neighbor].append(idr_neighbor)
                        fund_dom.append(idr_neighbor)
                        idrs_to_visit.appendleft(idr_neighbor)

    return FundDom(fund_dom, gens, edges_zipped)
