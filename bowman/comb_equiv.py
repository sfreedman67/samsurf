from collections import deque, namedtuple


class CombEquiv(namedtuple("CombEquiv", ["perm", "shift", "source", "target"])):
    __slots__ = ()

    def __repr__(self):
        reps = {(idx, 0): (self.perm[idx], self.shift[idx]) for idx in range(len(self.perm))}
        return f"CombEquiv({reps})"


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


# TODO: Move to triangulation?
def canonical_relabel(t, tri=0, edge=0):
    relabeling = {(tri, (edge + k) % 3): (0, k) for k in range(3)}
    tris_marked_to_visit = deque([(tri, edge)])
    tris_visited = {tri}
    while len(tris_visited) < len(t.triangles):
        tri, edge = tris_marked_to_visit.pop()
        for edge in [(edge + k) % 3 for k in range(3)]:
            tri_nbr, edge_nbr = t.gluings[(tri, edge)]
            if tri_nbr not in tris_visited:
                relabeling.update({
                    (tri_nbr, (edge_nbr + k) % 3): (len(tris_visited), (relabeling[(tri, edge)][1] + k) % 3)
                    for k in range(3)})
                tris_visited.add(tri_nbr)
                tris_marked_to_visit.appendleft((tri_nbr, (edge_nbr + 1) % 3))
    return relabeling


def generate_code_marked(t, tri, edge):
    relabeling = canonical_relabel(t, tri, edge)
    gluings_new = [(relabeling[e], relabeling[f])
                   for e, f in t.gluings.items()]
    gluings_new_ordered = [tuple(sorted(t)) for t in gluings_new]
    return hash(tuple(sorted(gluings_new_ordered))), (tri, edge)
