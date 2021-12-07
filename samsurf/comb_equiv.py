# TODO: Move to triangulation as a method?
def canonical_relabel(trin, tri_init, edge_init):
    """
    flag = (triangle, edge) pair
    INPUT:
    * trin - a triangulation
    """
    flags_to_visit = [(tri_init, edge_init)]
    relabelling = {}
    while flags_to_visit:
        (s, t) = flags_to_visit.pop()
        if (s, t) not in relabelling:
            relabelling[(s, t)] = divmod(len(relabelling), 3)
            neighbors = [trin.gluings[(s, t)], (s, (t + 2) % 3), (s, (t + 1) % 3)]
            flags_to_visit.extend(neighbors)
    return relabelling


def generate_code_marked(t, tri, edge):
    relabeling = canonical_relabel(t, tri, edge)
    gluings_new = [(relabeling[e], relabeling[f])
                   for e, f in t.gluings.items()]
    gluings_new_ordered = [tuple(sorted(t)) for t in gluings_new]
    return hash(tuple(sorted(gluings_new_ordered))), (tri, edge)
