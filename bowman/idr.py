from collections import namedtuple

import sage.all
from sage.all import *


class IDR(namedtuple("IDR", ["polygon", "labels_segment", "triangulation"])):
    __slots__ = ()

    def cross_segment(self, segment):
        hinges_degenerated = self.labels_segment[segment]

        triangulation_new = self.triangulation.flip_hinges(hinges_degenerated)

        return triangulation_new.IDR
