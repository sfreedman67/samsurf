from collections import namedtuple, defaultdict

import sage.all
from sage.all import *

from context import bowman
import bowman.halfplane
from bowman import halfplane


class IDR(namedtuple("IDR", ["polygon", "labels_segment", "triangulation"])):
    __slots__ = ()

    # TODO: should labels_segment be a list?

    def cross_segment(self, segment):
        hinges_degenerated = self.labels_segment[segment]

        triangulation_new = self.triangulation.flip_hinges(hinges_degenerated)

        return triangulation_new.IDR
