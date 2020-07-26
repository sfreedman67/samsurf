import sage.all
from sage.all import *

from context import bowman
import bowman.polygon
from bowman import polygon

def segments_real_intersect(s1, s2):
	_, (X1, Y1), (X2, Y2) = s1
	_, (X3, Y3), (X4, Y4) = s2

	I1 = (min(X1, X2), max(X1, X2))
	I2 = (min(X3, X4), max(X3, X4))

	Ia = (max(I1[0], I2[0]), min(I1[1], I2[1]))

	if max(X1, X2) < min(X3, X4):
		return False

	a1, b1, c1 = s1.halfplane
	a2, b2, c2 = s2.halfplane

	


def intersect_convex_polygons(P, Q):
	pass    

        