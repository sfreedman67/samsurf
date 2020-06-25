import bowman.halfplane
from bowman.halfplane import Edge

import bowman.trues_consecutive
from bowman.trues_consecutive import trues_consecutive


def truncate_ends(chain, point_head, point_tail):
    def process_edge(idx, edge):
        if idx == 0:
            return edge.chop_at(point_head)
        elif idx == len(chain) - 1:
            return edge.retract_to(point_tail)
        return edge

    return [process_edge(idx, edge) for idx, edge in enumerate(chain)]


def insert_new_halfplane(chain, halfplane):
    if not chain:
        start, end = halfplane.endpoints
        return [Edge(halfplane, start, end), Edge(None, end, start)]

    intersection_head = halfplane.edge_intersection(chain[-1])
    intersection_last = halfplane.edge_intersection(chain[0])
    edge_new = Edge(halfplane, intersection_head, intersection_last)

    return truncate_ends(chain, intersection_head, intersection_last) + [edge_new]


def intersect_halfplanes(halfplanes):
    if not halfplanes:
        return []

    edges = intersect_halfplanes(halfplanes[:-1])
    halfplane_curr = halfplanes[-1]

    intersects_curr = lambda edge: not halfplane_curr.is_edge_exterior(edge)
    edges_curr = list(trues_consecutive(edges, pred=intersects_curr))

    return insert_new_halfplane(edges_curr, halfplane_curr)
