# For each edge V1-V2 in the first polygon,
#     Let H := Half-plane tangenting V1-V2, with the remaining
#         vertices on the "inside".
#     Let C := New empty polygon.
#     For each edge V3-V4 in the second polygon,
#         Let X := The intersection between V3-V4 and H.
#         If V3 inside H, and V4 is outside H then,
#             Add V3 to C.
#             Add X to C.
#         Else if both V3 and V4 lies outside H then,
#             Skip.
#         Else if V3 outside H, and V4 is inside H then,
#             Add X to C.
#         Else
#             Add V3 to C.
#     Replace the second polygon with C.

import sage.all
from sage.all import *

from context import bowman
import bowman.polygon
from bowman import polygon

def intersect_convex_polygons(P1, P2):
    for edge in P1.edges:
        H = edge.halfplane
        C = Polygon()