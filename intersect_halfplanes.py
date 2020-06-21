import sage.all
from sage.all import *

import halfplane

def first_true(iterable, default=False, pred=None):
    """Returns the first true value in the iterable.

    If no true value is found, returns *default*

    If *pred* is not None, returns the first item
    for which pred(item) is true.

    """
    # first_true([a,b,c], x) --> a or b or c or x
    # first_true([a,b], x, f) --> a if f(a) else b if f(b) else x
    return next(filter(pred, iterable), default)

def intersect_halfplanes(halfplanes):
    if not halfplanes:
        return []
    elif len(halfplanes) == 1:
        assert False, "TODO: Decide a return type"

    intersection_current = intersect_halfplanes(halfplanes[:-1])

    # Test all the vertices:

    # Look for the preceeding edge

    # Look for the succeeding edge

    # retract pred
    # + NEW
    # + chop successor
    # + all edges from succ to pred