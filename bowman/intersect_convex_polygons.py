import itertools

import sage.all
from sage.all import *

from context import bowman
import bowman.polygon
from bowman import polygon
import bowman.triangulation
from bowman import triangulation


def segments_real_intersect(s1, s2):
    _, (X1, Y1), (X2, Y2) = s1
    _, (X3, Y3), (X4, Y4) = s2

    if max(X1, X2) < min(X3, X4):
        return None

    intersection = s1.halfplane.intersect_boundaries(s2.halfplane)

    if intersection is None:
        return None
    elif intersection.is_infinty:
        return polygon.Point(oo, QQ(0))
    elif intersection.u < max(min(X1, X2), min(X3, X4)):
        return None
    elif intersection.u > min(max(X1, X2), max(X3, X4)):
        return None
    return intersection


def intersect_convex_polygons(P, Q):
    edges_real_P = (edge for edge in P.edges if edge.halfplane is not None)
    for edge_P in edges_real_P:
        H = edge_P.halfplane

        edges_new = [component
                     for edge in Q.edges
                     for component in H.intersect_edge(edge)
                     if isinstance(component, polygon.Edge)]

        Q = polygon.Polygon._close_edge_chain(edges_new, H)

    return Q


def intersect_halfplanes(halfplanes):
    if not halfplanes:
        return polygon.Polygon([])
    elif len(halfplanes) == 1:
        [curr] = halfplanes
        return polygon.Polygon([polygon.Edge(curr, curr.start, curr.end),
                                polygon.Edge(None, curr.end, curr.start)])

    P = intersect_halfplanes(halfplanes[:len(halfplanes) // 2])
    Q = intersect_halfplanes(halfplanes[len(halfplanes) // 2:])

    if P is None or Q is None:
        return None

    return intersect_convex_polygons(P, Q)


if __name__ == "__main__":
    X = triangulation.Triangulation.arnoux_yoccoz(3)
    H = X.halfplanes

    P = intersect_halfplanes(H)
    P.plot().save("test.png")
