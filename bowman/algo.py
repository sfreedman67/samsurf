from collections import deque

from sage.all import *

from bowman import geom_equiv
from bowman.fund_dom import FundDom
from bowman.polygon import Point, Edge, Polygon
from bowman.idr import IDR


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
    if r0.has_self_equivalences:
        raise ValueError("Initial IDR has self-equivalences")
    edges_paired = {}
    edges_fund_dom = {x for x in r0.polygon}
    idrs_to_visit = deque([r0])
    fund_dom = {r0}
    generators = []
    code_to_idr = {r0.triangulation.code: r0}
    while idrs_to_visit:
        idr_curr = idrs_to_visit.pop()
        for idx_curr, edge_curr in enumerate(idr_curr.polygon):
            if edge_curr.reverse() not in edges_fund_dom and edge_curr not in edges_paired:
                trin_neighbor = idr_curr.get_trin_neighboring(idx_curr)
                code_neighbor = trin_neighbor.code
                idr_neighbor = trin_neighbor.idr

                if code_neighbor in code_to_idr:
                    ve, _ = geom_equiv.gen_geom_equiv(trin_neighbor, code_to_idr[code_neighbor].triangulation)
                    # first element of an equivalence is the matrix
                    m = sigma(ve)
                    edge_opp = edge_curr.apply_mobius(m).reverse()
                    edges_paired[edge_curr] = edge_opp
                    edges_paired[edge_opp] = edge_curr
                    generators.append(m)
                elif idr_neighbor.has_self_equivalences:
                    idr_neighbor_chopped, m_rot = chop_idr(idr_neighbor, edge_curr.reverse())
                    edge_in = idr_neighbor_chopped.polygon.edges[-2]
                    edge_out = idr_neighbor_chopped.polygon.edges[-1]

                    generators.append(sigma(m_rot))
                    edges_paired[edge_in] = edge_out
                    edges_paired[edge_out] = edge_in

                    code_to_idr[code_neighbor] = idr_neighbor_chopped
                    fund_dom.add(idr_neighbor_chopped)
                    edges_fund_dom |= {x for x in idr_neighbor_chopped.polygon}
                    idrs_to_visit.append(idr_neighbor_chopped)
                else:
                    code_to_idr[code_neighbor] = idr_neighbor
                    fund_dom.add(idr_neighbor)
                    edges_fund_dom |= {x for x in idr_neighbor.polygon}
                    idrs_to_visit.appendleft(idr_neighbor)

    if not _is_boundary_paired(fund_dom, edges_paired):
        raise BoundaryIsNotZippedPairwiseError

    return FundDom(fund_dom, generators, edges_paired)


def chop_idr(idr, edge):
    gen_rot = find_gen_ccw([m for m, _ in idr.triangulation.self_geom_equivs if m.trace().abs() < 2])
    center = fixed_point_elliptic(sigma(gen_rot))
    edge_translate = edge.apply_mobius(sigma(gen_rot))
    idx_edge = idr.polygon.edges.index(edge)
    idx_edge_translate = idr.polygon.edges.index(edge_translate)

    if idx_edge < idx_edge_translate:
        edges = idr.polygon.edges[idx_edge: idx_edge_translate]
    else:
        edges = idr.polygon.edges[idx_edge:] + idr.polygon.edges[: idx_edge_translate]

    edge_end_center = Edge.from_two_points(edge_translate.start, center)
    edge_center_start = Edge.from_two_points(center, edge.start)
    p_chopped = Polygon([*edges, edge_end_center, edge_center_start])
    idxs = [idr.polygon.edges.index(edge) for edge in edges]
    labels_new = {k: idr.labels_segment[idxs[k]] for k in range(len(edges))}
    idr_chopped = IDR(p_chopped, labels_new, idr.triangulation, folded=True)
    return idr_chopped, gen_rot


def find_gen_ccw(ms):
    e1 = sage.all.vector([1, 0])
    ms = [m for m in ms if sage.all.matrix([e1, m * e1]).determinant() > 0]
    return max(ms, key=order)


def order(m):
    if abs(m.trace()) >= 2:
        raise ValueError(f"{m} is not elliptic")
    k = 1
    A = m
    while A != sage.all.identity_matrix(2):
        A *= m
        k += 1
    return k
