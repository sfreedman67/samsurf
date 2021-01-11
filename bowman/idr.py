import sage.all
from sage.all import *

import collections

from bowman.algo import _find_veech_equivs


class IDR(collections.namedtuple("IDR", ["polygon", "labels_segment", "triangulation"])):
    __slots__ = ()

    def __repr__(self):
        return f"IDR with {len(self.polygon.edges)} sides"

    @property
    def is_trivial(self):
        return self.polygon is None or len(self.polygon.edges) <= 2

    def cross_segment(self, idx_segment):
        hinges_degenerated = self.labels_segment[idx_segment]

        triangulation_new = self.triangulation.flip_hinges(hinges_degenerated)

        return triangulation_new.idr

    @property
    def has_self_equivalences(self):
        return _find_veech_equivs(self, self) != [sage.all.identity_matrix(2)]

    @property
    def neighbors(self):
        return [self.cross_segment(k) for k in range(len(self.polygon.edges))]

    def plot(self):
        # TODO: add edge labels
        return self.polygon.plot()
