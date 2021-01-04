from collections import namedtuple, defaultdict

from sage.all import *


class FundDom(namedtuple("FundDom", ["idrs", "gens", "edges_zipped"])):
    __slots__ = ()

    def __repr__(self):
        return f"FundDom with {len(self.idrs)} idrs and {len(self.gens)} gens"

    @property
    def boundary(self):
        edges_total = {edge for idr in self.idrs for edge in idr.polygon.edges}
        boundary_unordered = {edge for edge in edges_total if edge.reverse() not in edges_total}
        return _order_boundary(boundary_unordered)

    @property
    def num_edges_folded(self):
        return sum(self.edges_zipped[edge] == edge.reverse()
                   for edge in self.boundary)

    @property
    def num_edges_normal(self):
        return QQ(1 / 2) * (len(self.boundary) - self.num_edges_folded)

    @property
    def cone_angles(self):
        matchings_edge = {edge: self.edges_zipped[edge] for edge in self.boundary}
        matchings_vertex = {e1.end: e2.end for e1, e2 in matchings_edge.items()}

        angles_boundary = {e1.end: e1.angle(e2)
                           for e1, e2 in zip(self.boundary,
                                             self.boundary[1:] + self.boundary[:1])}
        return [sum(angles_boundary[vertex] for vertex in orbit)
                for orbit in _get_orbits(matchings_vertex)]

    @property
    def area(self):
        return pi * (len(self.boundary) - 2) - sum(self.cone_angles)

    @property
    def num_cusps(self):
        return sum(bool(angle == 0) for angle in self.cone_angles)

    @property
    def points_orbifold(self):
        return [pi] * self.num_edges_folded + [angle for angle in self.cone_angles
                                               if RR(angle) != 0 and RR(angle / pi).round() != 2]

    @property
    def chi_top(self):
        return len(self.cone_angles) - self.num_edges_normal + 1

    @property
    def chi_orb(self):
        num_vertices = sum(a / (2 * pi) for a in self.cone_angles) + QQ(1 / 2) * self.num_edges_folded
        num_edges = self.num_edges_folded + self.num_edges_normal

        return num_vertices - num_edges + 1

    @property
    def genus(self):
        return (2 - self.chi_top) / 2

    @property
    def codes_to_idrs(self):
        codes_to_idrs = defaultdict(list)
        for idr in self.idrs:
            codes_to_idrs[idr.triangulation.code_comb].append(idr)
        return codes_to_idrs

    @property
    def proportions(self):
        return {code: RR(sum(idr.polygon.area for idr in idrs) / self.area)
                for code, idrs in self.codes_to_idrs.items()}

    def plot(self):
        return sum(r.plot() for r in self.idrs)


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
