from collections import deque, namedtuple


class CombEquiv(namedtuple("CombEquiv", ["perm", "shift", "source", "target"])):
    __slots__ = ()

    def __repr__(self):
        reps = {(idx, 0): (self.perm[idx], self.shift[idx]) for idx in range(len(self.perm))}
        return f"CombEquiv({reps})"

    @classmethod
    def from_edge(cls, t1, t2, e):
        matching = _develop_matching(t1, t2, e)
        if matching is not None:
            perm = tuple(matching[(k, 0)][0] for k in range(len(t1.triangles)))
            shifts = tuple(matching[(k, 0)][1] for k in range(len(t1.triangles)))
            return CombEquiv(perm, shifts, t1, t2)

    def move(self, edge):
        idx_tri, idx_edge = edge
        return self.perm[idx_tri], (idx_edge + self.shift[idx_tri]) % 3

    def respects_gluing(self, edge):
        image = self.move(edge)
        image_opp = self.target.gluings[image]
        edge_opp = self.source.gluings[edge]
        edge_opp_image = self.move(edge_opp)
        return image_opp == edge_opp_image

    @property
    def respects_gluings(self):
        """Check if e1~e2 in t1 implies F(e1) ~ F(e2) in t2"""
        return all(self.respects_gluing(edge) for edge in self.source.edges)


def _develop_matching(t1, t2, im00):
    queue = deque([(0, 0)])
    matching = {(0, 0): im00}

    while queue:
        tri, edge = queue.pop()
        edge_succ = (edge + 1) % 3
        edge_pred = (edge + 2) % 3
        tri_opp, edge_opp = t1.gluings[(tri, edge)]
        tri_im, edge_im = matching[(tri, edge)]

        for e1, e2 in [((tri, edge_succ), (tri_im, (edge_im + 1) % 3)),
                       ((tri, edge_pred), (tri_im, (edge_im + 2) % 3)),
                       ((tri_opp, edge_opp), t2.gluings[(tri_im, edge_im)])]:
            if e1 in matching:
                if matching[e1] != e2:
                    return None
            else:
                matching[e1] = e2
                queue.appendleft(e1)

    return matching


def gen_comb_equivs(t1, t2):
    poss_matchings = [CombEquiv.from_edge(t1, t2, (idx_tri, idx_edge))
                      for idx_tri in range(len(t1.triangles))
                      for idx_edge in range(3)]

    ces = [x for x in poss_matchings if x is not None]
    return ces


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
    return hash(tuple(sorted(gluings_new_ordered)))
