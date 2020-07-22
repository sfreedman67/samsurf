import itertools
from collections import defaultdict, namedtuple, deque

import sage.all
from sage.all import *
import flatsurf

from context import bowman
import bowman.halfplane
from bowman import halfplane


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

    def _constant_term(self):
        pass

    def _linear_term(self):
        pass

    def _quadratic_term(self):
        pass

    @property
    def halfplane(self):
        M_a = matrix([[x, y, y**2] for x, y in self.vectors])
        M_b = matrix([[x, y, x * y] for x, y in self.vectors])
        M_c = matrix([[x, y, x**2] for x, y in self.vectors])

        a = M_a.determinant()
        b = 2 * M_b.determinant()
        c = M_c.determinant()

        try:
            return halfplane.HalfPlane.from_ineq(a, b, c)
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
    def boundary_IDs(self):
        label_tri, label_edge = self.id_edge
        label_tri_opp, label_edge_opp = self.id_edge_opp

        return ((label_tri, (label_edge + 1) % 3),
                (label_tri, (label_edge + 2) % 3),
                (label_tri_opp, (label_edge_opp + 1) % 3),
                (label_tri_opp, (label_edge_opp + 2) % 3))


class Triangulation(namedtuple("Triangulation", ["triangles", "gluings", "field"])):

    @classmethod
    def _from_flatsurf(cls, X):
        # TODO: hardcode in the answers from sage so we can remove flatsurf

        DT = X.delaunay_triangulation()

        DT_polygons = [DT.polygon(i) for i in range(DT.num_polygons())]

        triangles = [Triangle(*[vector(edge) for edge in polygon.edges()])
                     for polygon in DT_polygons]

        gluings = {edge[0]: edge[1] for edge in DT.edge_iterator(gluings=True)}

        ring = DT.base_ring()

        return Triangulation(triangles, gluings, ring)

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
        halfplanes = [hinge.halfplane for hinge in self.hinges
                      if hinge.halfplane is not None]

        return list(set(halfplanes))

    @property
    def halfplane_to_id_edge(self):
        halfplane_to_id_edge = defaultdict(list)

        for hinge in self.hinges:
            halfplane = hinge.halfplane
            if halfplane is not None:
                halfplane_to_id_edge[halfplane].append(hinge.id_edge)

        return halfplane_to_id_edge

    def _triangles_after_flip(self, id_edge):
        triangles_new = self.triangles.copy()

        hinge_flipped = Hinge._from_id_edge(self, id_edge).flip()

        triangles_new[hinge_flipped.id_edge[0]] = hinge_flipped.triangle
        triangles_new[hinge_flipped.id_edge_opp[
            0]] = hinge_flipped.triangle_opp

        return triangles_new

    def _gluings_after_flip(self, id_edge):
        # TODO: fix this disaster

        IDs = Hinge._from_id_edge(self, id_edge).boundary_IDs
        new_name = {old_name: new_name
                    for new_name, old_name in zip(IDs,
                                                  IDs[1:] + IDs[:1])}

        substitution = lambda edge: edge if edge not in new_name.keys() else new_name[
            edge]

        return {substitution(key): substitution(value)
                for key, value in self.gluings.items()}

    def flip_hinge(self, id_edge):
        return Triangulation(self._triangles_after_flip(id_edge),
                             self._gluings_after_flip(id_edge),
                             self.field)

    def flip_hinges(self, ids_edges):
        triangulation = self
        for id_edge in ids_edges:
            triangulation = triangulation.flip_hinge(id_edge)
        return triangulation

    def plot_halfplanes(self, count=None):
        # TODO: Labels?
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
        return halfplane.HalfPlane.intersect_halfplanes(self.halfplanes)

    def iso_delaunay_complex(self, depth):
        IDRs_visited = set()
        IDRs_visited.add(self.IDR)

        queue = deque()
        queue.appendleft((self, self.IDR))

        # fig = self.IDR.plot()

        count = 0

        while(queue and count < depth):
            triangulation, IDR = queue[-1]
            for segment in IDR.edges:
                edges_degenerated = triangulation.halfplane_to_id_edge[
                    segment.halfplane]
                triangulation_new = triangulation.flip_hinges(
                    edges_degenerated)
                IDR_new = triangulation_new.IDR

                if IDR_new not in IDRs_visited:
                    # fig += IDR_new.plot()
                    count += 1
                    IDRs_visited.add(IDR_new)
                    queue.appendleft((triangulation_new, IDR_new))
            queue.pop()

        # fig.save("iso_delaunay_complex.png")
        # fig.show(xmin=-2, xmax=2, ymax=2)
