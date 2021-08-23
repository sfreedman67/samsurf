from collections import deque

from sage.all import *

from bowman import geom_equiv
from bowman.fund_dom import FundDom
from bowman.polygon import Point, Edge


def sigma(mat):
    [[a, b], [c, d]] = mat
    return sage.all.matrix([[a, -b], [-c, d]])


def fixed_point_elliptic(mat):
    [[a, b], [c, d]] = mat
    delta = mat.trace() ** 2 - 4 * mat.determinant()
    u = (a - d) / (2 * c)
    v2 = (delta / (4 * c ** 2)).abs()
    return Point(u, v2)


def _is_boundary_paired(idrs, pairings):
    edges_total = {x for idr in idrs for x in idr.polygon}
    boundary = {x for x in edges_total if x.reverse() not in edges_total}
    boundary_is_paired = all(x in pairings or x.reverse() in pairings for x in boundary)
    pairings_contains_boundary = all(x in boundary or x.reverse() in boundary for x in boundary)
    return boundary_is_paired and pairings_contains_boundary


class BoundaryIsNotZippedPairwiseError(Exception):
    pass


def generators_veech(trin):
    if not trin.is_delaunay:
        raise ValueError("Need to start with Delaunay Triangulation")
    r0 = trin.idr
    if r0.is_trivial:
        raise ValueError("The starting IDR is degenerate")

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
                trin_neighbor = idr_curr.get_trin_neighboring(idx)
                code_neighbor = trin_neighbor.code
                idr_neighbor = trin_neighbor.idr

                for edge in idr_neighbor.polygon:
                    if edge.start != oo and edge.start.u.C != 0:
                        if not edge.start.u.C.is_square():
                            assert False
                    if edge.end != oo and edge.end.u.C != 0:
                        if not edge.end.u.C.is_square():
                            assert False

                if code_neighbor in code_to_idr:
                    ve, _ = geom_equiv.gen_geom_equiv(trin_neighbor, code_to_idr[code_neighbor].triangulation)
                    # first element of an equivalence is the matrix
                    m = sigma(ve)
                    edge_opp = edge_curr.apply_mobius(m).reverse()
                    edges_paired[edge_curr] = edge_opp
                    edges_paired[edge_opp] = edge_curr
                    generators.append(m)
                elif idr_neighbor.has_self_equivalences:
                    m_rot = next(m for m, _ in idr_neighbor.triangulation.get_self_geom_equivs()
                                 if abs(m.trace()) < 2
                                 and edge_curr.end.apply_mobius(sigma(m)) == edge_curr.start)
                    e0 = edge_curr.reverse()
                    idr_center = fixed_point_elliptic(sigma(m_rot))
                    u0, u1 = edge_curr.start.u, idr_center.u
                    raise ValueError(f"IDR has self-equivalences")
                else:
                    code_to_idr[code_neighbor] = idr_neighbor
                    fund_dom.add(idr_neighbor)
                    edges_fund_dom |= {x for x in idr_neighbor.polygon}
                    idrs_to_visit.appendleft(idr_neighbor)

    if not _is_boundary_paired(fund_dom, edges_paired):
        raise BoundaryIsNotZippedPairwiseError

    return FundDom(fund_dom, generators, edges_paired)
