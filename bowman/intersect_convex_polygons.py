import sage.all
from sage.all import *

from context import bowman
import bowman.polygon
from bowman import polygon


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


def cross(v1, v2):
    return v1[0] * v2[1] - v2[0] * v1[1]


def intersect_convex_polygons(P, Q):
	count = 0
	p = P.edges[0]
	q = Q.edges[0]
	while(count < 2 * (len(P.edges) + len(Q.edges))):
		count += 1
		
	
	if Q.contains_point(p.start):
		return P
	elif P.contains_point(q.start):
		return Q
	return None

