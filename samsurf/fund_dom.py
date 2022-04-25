from collections import defaultdict, Counter

from sage.all import *


class FundDom:
    def __init__(self, idrs, gens, pairings):
        self.idrs = idrs
        self.gens = gens
        self.pairings = pairings

        self._boundary = None

    def __repr__(self):
        return f"FundDom with {len(self.idrs)} idrs and {len(self.gens)} gens"

    def __len__(self):
        return len(self.idrs)

    @property
    def boundary(self):
        if self._boundary is not None:
            return self._boundary
        edges_total = {edge for idr in self.idrs for edge in idr.polygon.edges}
        boundary_unordered = {edge for edge in edges_total if edge.reverse() not in edges_total}
        self._boundary = _order_boundary(boundary_unordered)
        return self._boundary

    @property
    def num_edges_folded(self):
        return sum(self.pairings[edge] == edge
                   for edge in self.boundary)

    @property
    def num_edges_normal(self):
        return QQ(1 / 2) * (len(self.boundary) - self.num_edges_folded)

    @property
    def cone_angles(self):
        matchings_vertex = {edge.start: self.pairings[edge].end for edge in self.boundary}
        angles_boundary = {e1.end: e1.angle(e2)
                           for e1, e2 in zip(self.boundary,
                                             self.boundary[1:] + self.boundary[:1])}
        return [sum(angles_boundary[vertex] for vertex in orbit)
                for orbit in _get_orbits(matchings_vertex)]

    @property
    def cusps(self):
        return sum(bool(angle == 0) for angle in self.cone_angles)

    @property
    def points_orbifold(self):
        orders = [RR(RR(2 * pi) / angle).round() for angle in self.cone_angles if angle != RR(0)]
        return [2] * self.num_edges_folded + [x for x in orders if x != 1]

    @property
    def chi_top(self):
        return len(self.cone_angles) - self.num_edges_normal + 1

    @property
    def chi_orb(self):
        num_vertices = sum(a for a in self.cone_angles) / (2 * RR(pi)) + QQ(1 / 2) * self.num_edges_folded
        num_edges = self.num_edges_folded + self.num_edges_normal

        return RR(num_vertices - num_edges + 1).nearby_rational(max_error=0.00001)

    @property
    def genus(self):
        return (2 - self.chi_top) / 2

    @property
    def codes_comb_to_idrs(self):
        codes_to_idrs = defaultdict(list)
        for idr in self.idrs:
            codes_to_idrs[idr.triangulation.code_comb].append(idr)
        return codes_to_idrs

    @property
    def proportions(self):
        areas = Counter()
        for idr in self.idrs:
            areas[idr.triangulation.code_comb] += (idr.area / self.area)
        return areas

    def plot(self):
        return sum(x.plot() for x in self.idrs)

    @property
    def area(self):
        return sum(p.area for p in self.idrs)


def _order_boundary(boundary):
    edge1 = next(iter(boundary))
    boundary_ordered = [edge1]
    while len(boundary_ordered) < len(boundary):
        edge2 = next(edge for edge in boundary if edge.start == edge1.end)
        boundary_ordered.append(edge2)
        edge1 = edge2
    return boundary_ordered


def _get_orbits(f):
    seen = set()
    orbits = []
    for p in f:
        if p not in seen:
            orbit = _get_orbit(f, p)
            seen |= set(orbit)
            orbits.append(orbit)
    return orbits


def _get_orbit(f, p):
    orbit = [p]
    fp = f[p]
    while fp not in orbit:
        orbit.append(fp)
        fp = f[fp]
    return orbit
