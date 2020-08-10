import itertools
from collections import defaultdict, namedtuple, deque
from functools import lru_cache

import sage.all
from sage.all import *
import flatsurf

from context import bowman
import bowman.halfplane
from bowman import halfplane
import bowman.idr
from bowman import idr


class Triangle(namedtuple("Triangle", ["v0", "v1", "v2"])):
    __slots__ = ()

    def __new__(cls, v0, v1, v2):
        if sum((v0, v1, v2)) != vector([0, 0]):
            raise ValueError("sides do not close up")
        elif matrix([v0, -v2]).determinant() <= 0:
            raise ValueError("sides are not oriented correctly")

        self = super(Triangle, cls).__new__(cls, v0, v1, v2)
        return self

    def apply_matrix(self, M):
        return Triangle(*(M * side for side in self))


class Hinge(namedtuple("Hinge", ["vectors", "id_edge", "id_edge_opp"])):
    __slots__ = ()

    @property
    def coords(self):
        return tuple(coord for vector in self.vectors for coord in vector)

    def __hash__(self):
        return hash((self.coords,))

    @classmethod
    def _from_id_edge(cls, trin, id_edge):
        label_tri, label_edge = id_edge

        id_edge_opp = trin.gluings[id_edge]

        label_tri_opp, label_edge_opp = id_edge_opp

        tri = trin.triangles[label_tri]
        tri_opp = trin.triangles[label_tri_opp]

        edge = tri[label_edge]
        edge_opp = tri_opp[label_edge_opp]

        if edge != -edge_opp:
            raise ValueError("Edges either nonparallel or improperly oriented")

        v0 = tri[(label_edge + 1) % 3]
        v1 = edge_opp
        v2 = -tri_opp[(label_edge_opp - 1) % 3]

        return Hinge((v0, v1, v2), id_edge, id_edge_opp)

    def flip(self):
        v0, v1, v2 = self.vectors
        return Hinge((v1 - v0, v2 - v0, -v0), self.id_edge, self.id_edge_opp)

    @property
    def incircle_det(self):
        """(p2 is inside/on/outisde oriented circle 0-P0-P1) iff (det </==/> 0) """
        return matrix([[x, y, x**2 + y**2] for x, y in self.vectors]).determinant()

    @property
    @lru_cache(None)
    def _coefficients(self):
        (x0, y0), (x1, y1), (x2, y2) = self.vectors

        m02 = x1 * y2 - x2 * y1
        m12 = x0 * y2 - x2 * y0
        m22 = x0 * y1 - x1 * y0

        a = y0**2 * m02 - y1**2 * m12 + y2**2 * m22
        b = 2 * (x0 * y0 * m02 - x1 * y1 * m12 + x2 * y2 * m22)
        c = x0**2 * m02 - x1**2 * m12 + x2**2 * m22

        return (a, b, c)

    @property
    def halfplane(self):
        try:
            return halfplane.HalfPlane.from_ineq(*self._coefficients)
        except ValueError:
            return None

    @property
    def triangle(self):
        v0, v1, v2 = self.vectors

        sides_ordered = sorted([(self.id_edge[1], -v1),
                                ((self.id_edge[1] + 1) % 3, v0),
                                ((self.id_edge[1] + 2) % 3, v1 - v0)])
        return Triangle(*(vector for _, vector in sides_ordered))

    @property
    def triangle_opp(self):
        v0, v1, v2 = self.vectors

        sides_ordered = sorted([(self.id_edge_opp[1], v1),
                                ((self.id_edge_opp[1] + 1) % 3, v2 - v1),
                                ((self.id_edge_opp[1] + 2) % 3, -v2)])

        return Triangle(*(vector for _, vector in sides_ordered))

    @property
    def _IDs_boundary(self):
        '''return the edge IDs of the boundary of the hinge
        starting in the NE and moving Clockwise'''

        label_tri, label_edge = self.id_edge
        label_tri_opp, label_edge_opp = self.id_edge_opp

        SE = (label_tri, (label_edge + 1) % 3)
        NE = (label_tri, (label_edge + 2) % 3)
        NW = (label_tri_opp, (label_edge_opp + 1) % 3)
        SW = (label_tri_opp, (label_edge_opp + 2) % 3)

        return (NE, SE, SW, NW)


class Triangulation(namedtuple("Triangulation", ["triangles", "gluings", "field"])):

    @classmethod
    def _from_flatsurf(cls, X):
        DT = X.delaunay_triangulation()

        DT_polygons = [DT.polygon(i) for i in range(DT.num_polygons())]

        triangles = [Triangle(*[vector(edge) for edge in polygon.edges()])
                     for polygon in DT_polygons]

        gluings = {edge[0]: edge[1] for edge in DT.edge_iterator(gluings=True)}

        field = DT.base_ring()

        return Triangulation(triangles, gluings, field)

    @classmethod
    def square_torus(cls):
        return cls._from_flatsurf(flatsurf.translation_surfaces.square_torus())

    @classmethod
    def regular_octagon(cls):
        return cls._from_flatsurf(flatsurf.translation_surfaces.regular_octagon())

    @classmethod
    def arnoux_yoccoz(cls, g):
        if g < 3:
            raise ValueError("g must be >= 3")
        return cls._from_flatsurf(flatsurf.translation_surfaces.arnoux_yoccoz(g))

    @classmethod
    def octagon_and_squares(cls):
        return cls._from_flatsurf(flatsurf.translation_surfaces.octagon_and_squares())

    @property
    def edges(self):
        num_triangles = len(self.triangles)
        edges = itertools.product(range(num_triangles), range(3))
        return [edge for edge in edges if edge < self.gluings[edge]]

    def apply_matrix(self, M):
        return Triangulation([triangle.apply_matrix(M) for triangle in self.triangles], self.gluings, self.field)

    @property
    def hinges(self):
        return [Hinge._from_id_edge(self, edge) for edge in self.edges]

    def is_delaunay(self, non_degenerate=False):
        return all(hinge.incircle_det > 0
                   if non_degenerate
                   else hinge.incircle_det >= 0
                   for hinge in self.hinges)

    @property
    def halfplanes(self):
        return list(filter(lambda x: x is not None,
                           (hinge.halfplane for hinge in self.hinges)))

    @property
    def _halfplanes_to_hinges_degen(self):
        halfplanes_labelled = [(hinge.halfplane, hinge.id_edge)
                               for hinge in self.hinges]

        dd = defaultdict(list)
        for halfplane, id_hinge in halfplanes_labelled:
            dd[halfplane].append(id_hinge)

        # remove degen halfplanes that are None
        dd.pop(None, None)

        return dd

    def _triangles_after_flip(self, hinge_flipped):
        idx_tri_new = hinge_flipped.id_edge[0]
        idx_tri_opp_new = hinge_flipped.id_edge_opp[0]

        triangles_new = []
        for idx, triangle in enumerate(self.triangles):
            if idx == idx_tri_new:
                triangles_new.append(hinge_flipped.triangle)
            elif idx == idx_tri_opp_new:
                triangles_new.append(hinge_flipped.triangle_opp)
            else:
                triangles_new.append(triangle)

        return triangles_new

    def _id_edge_after_flip(self, hinge_flipped, edge):
        NE, SE, SW, NW = hinge_flipped._IDs_boundary

        IDs_new = {NE: SE, SE: SW, SW: NW, NW: NE}

        if edge in IDs_new:
            return IDs_new[edge]
        return edge

    def _gluings_after_flip(self, hinge_flipped):
        return {self._id_edge_after_flip(hinge_flipped, key):
                self._id_edge_after_flip(hinge_flipped, value)
                for key, value in self.gluings.items()}

    def flip_hinge(self, id_edge):
        hinge_flipped = Hinge._from_id_edge(self, id_edge).flip()

        return Triangulation(self._triangles_after_flip(hinge_flipped),
                             self._gluings_after_flip(hinge_flipped),
                             self.field)

    def flip_hinges(self, ids_edges):
        triangulation = self
        for id_edge in ids_edges:
            triangulation = triangulation.flip_hinge(id_edge)
        return triangulation

    def plot_halfplanes(self, count=None):
        figure = sum(itertools.islice((halfplane.plot()
                                       for halfplane in self.halfplanes), count))
        if count is not None:
            plt_final = figure[-1]
            opt = plt_final.options()
            opt["linestyle"] = "--"
            plt_final.set_options(opt)

        return figure

    @property
    def IDR(self):
        halfplane_to_ids_hinge = self._halfplanes_to_hinges_degen
        halfplanes = list(halfplane_to_ids_hinge.keys())

        P = halfplane.HalfPlane.intersect_halfplanes(halfplanes)

        labels_segment = {idx: halfplane_to_ids_hinge[segment.halfplane]
                          for idx, segment in enumerate(P.edges)}

        return idr.IDR(P, labels_segment, self)

    def iso_delaunay_complex(self, limit):
        IDR_start = self.IDR

        polygons_visited = {IDR_start.polygon}
        segments_crossed = set()

        queue = deque([IDR_start])

        while len(polygons_visited) < limit:
            IDR = queue.pop()
            segments_uncrossed = [(idx, segment)
                                  for idx, segment in enumerate(IDR.polygon.edges)
                                  if not segment in segments_crossed]

            for idx_segment, segment in segments_uncrossed:
                IDR_new = IDR.cross_segment(idx_segment)
                segments_crossed |= {segment, segment.reverse()}

                if IDR_new.polygon not in polygons_visited:
                    polygons_visited.add(IDR_new.polygon)
                    queue.appendleft(IDR_new)

        return list(polygons_visited)
