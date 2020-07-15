from context import bowman

import bowman.polygon
from bowman.polygon import Edge


def intersect_halfplanes(halfplanes):
    if not halfplanes:
        return []

    polygon_previous = intersect_halfplanes(halfplanes[:-1])
    current = halfplanes[-1]

    if polygon_previous is None:
        return None

    elif polygon_previous == []:
        return [Edge(current, current.start, current.end),
                Edge(None, current.end, current.start)]

    return bowman.polygon.intersect_polygon_halfplane(polygon_previous, current)
