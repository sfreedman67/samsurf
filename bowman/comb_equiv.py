from collections import deque, namedtuple


class CombEquiv(namedtuple("CombEquiv", ["perm", "shift", "source", "target"])):
    __slots__ = ()

    def __repr__(self):
        return f"CombEquiv(perm={self.perm}, shift={self.shift})"

    @classmethod
    def _from_matching_partial(cls, matching_partial, source, target):
        perm = []
        shifts = []

        for edge in [(k, 0) for k in range(len(source.triangles))]:
            lab_tri_im, lab_edge_im = matching_partial[edge]
            perm.append(lab_tri_im)
            shifts.append((lab_edge_im - edge[1]) % 3)

        return CombEquiv(tuple(perm), tuple(shifts), source, target)

    @classmethod
    def from_edge(cls, t1, t2, e):
        matching_partial = _get_matching_partial_from_edge(t1, t2, e)

        return CombEquiv._from_matching_partial(*matching_partial)

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


def _get_matching_partial_from_edge(t1, t2, edge_00_im):
    def extend_to_nbrs(e, f):
        return {(e[0], (e[1] + k) % 3): (f[0], (f[1] + k) % 3)
                for k in range(3)}

    tris_to_visit = deque([0])
    matching_partial = extend_to_nbrs((0, 0), edge_00_im)

    tris_visited = {0}
    while len(tris_visited) != len(t1.triangles):
        tri_curr = tris_to_visit.pop()
        edges_curr = [(tri_curr, k) for k in range(3)]
        for edge in edges_curr:
            edge_nbr = t1.gluings[edge]
            edge_im = matching_partial[edge]
            edge_im_nbr = t2.gluings[edge_im]
            tri_nbr = edge_nbr[0]
            if tri_nbr not in tris_visited:
                matching_partial.update(extend_to_nbrs(edge_nbr, edge_im_nbr))
                tris_visited.add(tri_nbr)
                tris_to_visit.appendleft(tri_nbr)

    return matching_partial, t1, t2


def gen_comb_equivs(t1, t2):
    poss_matchings = [CombEquiv.from_edge(t1, t2, (idx_tri, idx_edge))
                      for idx_tri in range(len(t1.triangles))
                      for idx_edge in range(3)]

    return [x for x in poss_matchings if x.respects_gluings]
