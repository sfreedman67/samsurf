import collections


class IDR(collections.namedtuple("IDR", ["polygon", "labels_segment", "triangulation"])):
    __slots__ = ()

    @property
    def is_nontrivial(self):
        return self.polygon is not None and len(self.polygon.edges) > 2

    def cross_segment(self, idx_segment):
        hinges_degenerated = self.labels_segment[idx_segment]

        triangulation_new = self.triangulation.flip_hinges(hinges_degenerated)

        return triangulation_new.get_idr

    @property
    def neighbors(self):
        return [self.cross_segment(k) for k in range(len(self.polygon.edges))]

    def plot(self):
        return self.polygon.plot()
