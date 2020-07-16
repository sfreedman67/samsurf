import itertools
from collections import defaultdict, namedtuple

import sage.all
from sage.all import *
import flatsurf

from context import bowman
import bowman.halfplane
from bowman import halfplane


class Triangle:
    # TODO: Named Tuple

    def __init__(self, edges):
        if sum(edges) != vector([0, 0]):
            raise ValueError("sides do not close up")
        elif matrix([edges[0], -edges[2]]).determinant() <= 0:
            raise ValueError("sides are not oriented correctly")
        else:
            self.edges = edges

    def __repr__(self):
        return f'Triangle{self.edges}'

    def __eq__(self, other):
        if isinstance(other, Triangle):
            return self.edges == other.edges
        return NotImplemented

    def apply_matrix(self, M):
        return Triangle([M * edge for edge in self.edges])


# class Hinge(namedtuple("Hinge", ["v0, v1, v2, id"])):
class Hinge:
    # TODO: Named Tuple

    def __init__(self, v0, v1, v2):
        self.vectors = (v0, v1, v2)

    def __repr__(self):
        return f"Hinge{self.vectors}"

    def __eq__(self, other):
        if isinstance(other, Hinge):
            return self.vectors == other.vectors
        return NotImplemented

    @classmethod
    def _from_triangles(cls, triangle1, label_e1, triangle2, label_e2):
        edge1, edge2 = triangle1.edges[label_e1], triangle2.edges[label_e2]

        if edge1 != -edge2:
            raise ValueError(
                "Edges are either nonparallel or oriented incorrectly")

        v0 = triangle2.edges[(label_e2 + 1) % 3]
        v1 = edge1
        v2 = -triangle1.edges[(label_e1 - 1) % 3]

        return Hinge(v0, v1, v2)

    @property
    def incircle_det(self):
        """(p2 is inside/on/outisde oriented circle 0-P0-P1) iff (det </==/> 0) """
        return matrix([[x, y, x**2 + y**2] for x, y in self.vectors]).determinant()

    @property
    def halfplane(self):
        M_a = matrix([[x, y, y**2] for x, y in self.vectors])
        M_b = matrix([[x, y, x * y] for x, y in self.vectors])
        M_c = matrix([[x, y, x**2] for x, y in self.vectors])

        a = M_a.determinant()
        b = 2 * M_b.determinant()
        c = M_c.determinant()

        try:
            H = halfplane.HalfPlane.from_ineq(a, b, c)
        except ValueError:
            H = None

        return H


class Triangulation:
    # TODO: Named tuple? Assuming I never change gluings.
    # TODO: Use a tuple of triangles to prevent side effects

    def __init__(self, triangles, gluings, base_ring):
        self.triangles = triangles
        self.gluings = gluings
        self.base_ring = base_ring

    def __repr__(self):
        return f'Triangulation(edges={self.triangles}, gluings={self.gluings})'

    def __eq__(self, other):
        if isinstance(other, Triangulation):
            return self.triangles == other.triangles and self.gluings == other.gluings
        return NotImplemented

    @classmethod
    def _from_flatsurf(cls, X):
        # TODO: hardcode in the answers from sage so we can remove flatsurf

        DT = X.delaunay_triangulation()

        DT_polygons = [DT.polygon(i) for i in range(DT.num_polygons())]

        triangles = [Triangle([vector(edge) for edge in polygon.edges()])
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

    def edges(self, gluings=False):
        num_triangles = len(self.triangles)
        edges = itertools.product(range(num_triangles), range(3))
        reps = [edge for edge in edges if edge < self.gluings[edge]]
        return reps if not gluings else [(rep, self.gluings[rep]) for rep in reps]

    def apply_matrix(self, M):
        return Triangulation([triangle.apply_matrix(M) for triangle in self.triangles], self.gluings, self.base_ring)

    def _hinge(self, edge):
        edge_opposite = self.gluings[edge]
        tri1 = self.triangles[edge[0]]
        tri2 = self.triangles[edge_opposite[0]]
        return Hinge._from_triangles(tri1, edge[1], tri2, edge_opposite[1])

    @property
    def hinges(self):
        return [self._hinge(edge) for edge in self.edges()]

    def is_delaunay(self, strict=True):
        has_valid_sign = lambda x: bool(x > 0) or (not strict and bool(x == 0))
        return all(has_valid_sign(hinge.incircle_det) for hinge in self.hinges)

    @property
    def halfplanes(self):
        halfplanes = [hinge.halfplane for hinge in self.hinges
                      if hinge.halfplane is not None]

        return list(set(halfplanes))

    @property
    def halfplanes_to_hinges(self):
        halfplanes_to_hinges = defaultdict(list)
        
        for hinge in self.hinges:
            halfplane = hinge.halfplane
            if halfplane is not None:
                halfplanes_to_hinges[halfplane].append(hinge)
        
        return halfplanes_to_hinges

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
