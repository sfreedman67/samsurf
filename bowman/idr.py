from sage.all import *

import collections

from bowman.geom_equiv import gen_geom_equivs


class IDR(collections.namedtuple("IDR", ["polygon", "labels_segment", "triangulation"])):
    __slots__ = ()

    def __repr__(self):
        return f"IDR with {len(self.polygon.edges)} sides"

    def __hash__(self):
        return hash(self.polygon)

    @property
    def is_trivial(self):
        return self.polygon is None or len(self.polygon.edges) <= 2

    def get_idr_neighboring(self, idx_side):
        return self.get_trin_neighboring(idx_side).idr

    def get_trin_neighboring(self, idx_side):
        hinges_degenerated = self.labels_segment[idx_side]
        return self.triangulation.flip_hinges(hinges_degenerated)

    @property
    def has_self_equivalences(self):
        """
        Checks if there is an equivalence from triangulation to itself
        whose matrix is not the identity or the 180 degree rotation
        #TODO: Whether this is the correct notion of self equivalence is a
        conversation for another time
        """
        self_geom_equiv_matrices = [equiv[0] for equiv in
                                    gen_geom_equivs(self.triangulation,
                                                    self.triangulation)]
        # first element of an equivalence is the matrix
        allowed_equiv_matrices = [matrix([[1, 0], [0, 1]]),
                                  matrix([[-1, 0], [0, -1]])]
        return any(x not in allowed_equiv_matrices
                   for x in self_geom_equiv_matrices)

    @property
    def neighbors(self):
        return [self.get_idr_neighboring(k) for k in range(len(self.polygon.edges))]

    @property
    def area(self):
        return self.polygon.area

    def plot(self):
        return self.polygon.plot()

    def contains_point(self, point):
        return self.polygon.contains_point(point)

    @property
    def cusps(self):
        cusps_list = []
        for edge in self.polygon.edges:
            if edge.start == oo:
                cusps_list.append(oo)
            elif sign(edge.start.v2) == 0:
                cusps_list.append(edge.start.u)
        return cusps_list
