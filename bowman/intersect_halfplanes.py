import bowman.polygon as polygon


def intersect_halfplanes(halfplanes):
    if not halfplanes:
        return []

    polygon_previous = intersect_halfplanes(halfplanes[:-1])
    current = halfplanes[-1]

    if polygon_previous is None:
        return None

    elif polygon_previous == []:
        return [polygon.Edge(current, current.start, current.end),
                polygon.Edge(None, current.end, current.start)]

    return polygon.intersect_polygon_halfplane(polygon_previous, current)

    
