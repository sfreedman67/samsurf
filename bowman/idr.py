from collections import namedtuple

import sage.all
from sage.all import *


class IDR(namedtuple("IDR", ["polygon", "labels_hinge", "triangulation"])):
    __slots__ = ()

    # TODO: access segments quicker?

    def cross_segment(self, segment):
        hinges_degenerated = self.labels_hinge[segment.halfplane]

        triangulation_new = self.triangulation.flip_hinges(hinges_degenerated)

        return triangulation_new.IDR
