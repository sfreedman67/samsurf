from collections import namedtuple


class CombEquiv(namedtuple("CombEquiv", ["perm", "shift", "source", "target"])):
    __slots__ = ()

    def __repr__(self):
        reps = {(idx, 0): (self.perm[idx], self.shift[idx]) for idx in range(len(self.perm))}
        return f"CombEquiv({reps})"


# TODO: Wow, do I really not need this anymore?
def gen_comb_equivs(t1, t2):
    if t1.code_comb != t2.code_comb:
        return []

    ces = []
    _, (tri1, edge1) = next(iter(t1.codes_comb))
    d1 = canonical_relabel(t1, tri1, edge1)
    for (tri2, edge2) in [(tri2, edge2) for _, (tri2, edge2) in t2.codes_comb]:
        d2 = canonical_relabel(t2, tri2, edge2)
        d2inv = {v: k for k, v in d2.items()}
        equiv = {(k, 0): d2inv[d1[(k, 0)]] for k in range(len(t1.triangles))}
        perm, shift = zip(*tuple(equiv[(k, 0)] for k in range(len(t2.triangles))))
        ces.append(CombEquiv(perm, shift, t1, t2))
    return ces


# TODO: Move to triangulation as a method?
def canonical_relabel(trin, tri_init, edge_init):
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
