from collections import namedtuple, defaultdict

import sage.all
from sage.all import *

from context import bowman
import bowman.halfplane
from bowman import halfplane


class IDR(namedtuple("IDR", ["polygon", "labels_segment", "triangulation"])):
    __slots__ = ()

    def cross_segment(self, idx_segment):
        hinges_degenerated = self.labels_segment[idx_segment]

        triangulation_new = self.triangulation.flip_hinges(hinges_degenerated)

        return triangulation_new.IDR
