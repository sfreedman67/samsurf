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

    def cross_segment(self, idx_segment):
        hinges_degenerated = self.labels_segment[idx_segment]

        triangulation_new = self.triangulation.flip_hinges(hinges_degenerated)

        return triangulation_new.idr

    @property
    def has_self_equivalences(self):
        return len({tuple(tuple(row) for row in x)
                    for x in gen_geom_equivs(self.triangulation, self.triangulation)}) > 1

    @property
    def neighbors(self):
        return [self.cross_segment(k) for k in range(len(self.polygon.edges))]

    @property
    def area(self):
        return self.polygon.area

    def plot(self):
        return self.polygon.plot()

    def contains_point(self, point):
        return self.polygon.contains_point(point)
