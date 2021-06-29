from collections import deque

from sage.all import *

from bowman import geom_equiv
from bowman.fund_dom import FundDom


def sigma(mat):
    [[a, b], [c, d]] = mat
    return sage.all.matrix([[a, -b], [-c, d]])


def _is_boundary_paired(idrs, pairings):
    edges_total = {x for idr in idrs for x in idr.polygon}
    boundary = {x for x in edges_total if x.reverse() not in edges_total}
    return all(pairings[x] in boundary for x in boundary)


class BoundaryIsNotZippedPairwiseError(Exception):
    pass


def generators_veech(trin):
    if not trin.is_delaunay:
        raise ValueError("Need to start with Delaunay Triangulation")
    r0 = trin.idr
    if r0.is_trivial:
        raise ValueError("The starting IDR is degenerate")

    elif r0.has_self_equivalences:
        raise ValueError("Starting IDR has self-equivalences")

    edges_paired = {}
    edges_fund_dom = {x for x in r0.polygon}
    idrs_to_visit = deque([r0])
    fund_dom = {r0}
    generators = []
    code_to_idr = {r0.triangulation.code: r0}
    while idrs_to_visit:
        idr_curr = idrs_to_visit.pop()
        for idx, edge_curr in enumerate(idr_curr.polygon):
            if edge_curr.reverse() not in edges_fund_dom and edge_curr not in edges_paired:
                idr_neighbor = idr_curr.cross_segment(idx)
                trin_neighbor = idr_neighbor.triangulation
                code_neighbor = trin_neighbor.code
                if code_neighbor in code_to_idr:
                    ve = geom_equiv.gen_geom_equiv(trin_neighbor, code_to_idr[code_neighbor].triangulation)
                    m = sigma(ve)
                    edge_opp = edge_curr.apply_mobius(m).reverse()
                    edges_paired[edge_curr] = edge_opp
                    edges_paired[edge_opp] = edge_curr
                    generators.append(m)
                else:
                    code_to_idr[code_neighbor] = idr_neighbor
                    fund_dom.add(idr_neighbor)
                    edges_fund_dom |= {x for x in idr_neighbor.polygon}
                    idrs_to_visit.appendleft(idr_neighbor)

    if not _is_boundary_paired(fund_dom, edges_paired):
        raise BoundaryIsNotZippedPairwiseError

    return FundDom(fund_dom, generators, edges_paired)
