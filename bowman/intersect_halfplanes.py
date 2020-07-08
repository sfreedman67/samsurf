import itertools
from itertools import takewhile, dropwhile, chain

from bowman.polygon import Edge


def intersect_halfplanes(halfplanes):
    if not halfplanes:
        return []

    intersection_previous = intersect_halfplanes(halfplanes[1:])
    current = halfplanes[0]

    if intersection_previous is None:
        return None

    elif intersection_previous == []:
        return [Edge(current, current.start, current.end),
                Edge(None, current.end, current.start)]

    curr_intersect_prev = [component for edge in intersection_previous
                           for component in current.intersect_edge(edge)
                           if isinstance(component, Edge)]

    if curr_intersect_prev == intersection_previous:
        return intersection_previous
    elif not curr_intersect_prev:
        return None


    [head_idx] = [idx for idx, edge in enumerate(curr_intersect_prev)
                  if current.contains_point(edge.start, on_boundary=True)]

    edge_chain = [*curr_intersect_prev[head_idx:],
                  *curr_intersect_prev[:head_idx]]

    edge_new = Edge(current, edge_chain[-1].end, edge_chain[0].start)

    return [*edge_chain, edge_new]
