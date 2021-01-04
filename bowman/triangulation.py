import itertools
from collections import defaultdict, namedtuple, deque
from functools import lru_cache

from sage.all import *
import flatsurf

from bowman import halfplane
from bowman import idr
from bowman import comb_equiv
from bowman import algo


class Triangle(namedtuple("Triangle", ["v0", "v1", "v2"])):
    __slots__ = ()

    def __new__(cls, v0, v1, v2):
        if sum((v0, v1, v2)) != 0:
            raise ValueError("sides do not close up")
        elif sage.all.matrix([v0, -v2]).determinant() <= 0:
            raise ValueError("sides are not oriented correctly")

        self = super(Triangle, cls).__new__(cls, v0, v1, v2)
        return self

    def reflect(self, idx):
        def reflect_vector(v, w):
            w_parallel = (v.dot_product(w) / v.dot_product(v)) * v
            w_perp = w - w_parallel
            return w - 2 * w_perp

        v_axis = self[idx]
        v_succ = self[(idx + 1) % 3]
        v_pred = self[(idx + 2) % 3]
        sides_new = {idx: -v_axis,
                     (idx + 1) % 3: reflect_vector(v_axis, -v_pred),
                     (idx + 2) % 3: reflect_vector(v_axis, -v_succ)}
        return Triangle(sides_new[0], sides_new[1], sides_new[2])

    def plot(self, basepoint=sage.all.zero_vector(2)):
        return sage.all.polygon2d(self.vertices(basepoint)).plot()

    def vertices(self, basepoint=sage.all.zero_vector(2)):
        return [basepoint, basepoint + self.v0, basepoint - self.v2]

    def __hash__(self):
        return hash(tuple(coord for vertex in self.vertices() for coord in vertex))


class Hinge(namedtuple("Hinge", ["vectors", "id_edge", "id_edge_opp"])):
    __slots__ = ()

    @property
    def coordinates(self):
        return tuple(coord for vector in self.vectors for coord in vector)

    def __hash__(self):
        return hash((self.coordinates,))

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

    @property
    def is_convex(self):
        v0, v1, v2 = self.vectors
        boundary = [v0, v1 - v0, v2 - v1, -v2]
        crosses = [w0x * w1y - w1x * w0y
                   for (w0x, w0y), (w1x, w1y)
                   in zip(boundary, boundary[1:] + boundary[:1])]

        all_positive = all(bool(cross > 0) for cross in crosses)
        all_negative = all(bool(cross < 0) for cross in crosses)

        return all_positive or all_negative

    def flip(self):
        v0, v1, v2 = self.vectors
        return Hinge((v1 - v0, v2 - v0, -v0), self.id_edge, self.id_edge_opp)

    @property
    def incircle_det(self):
        """(p2 is inside/on/outisde oriented circle 0-P0-P1) iff (det </==/> 0) """
        return sage.all.matrix([[x, y, x ** 2 + y ** 2] for x, y in self.vectors]).determinant()

    @property
    @lru_cache(None)
    def _coefficients(self):
        (x0, y0), (x1, y1), (x2, y2) = self.vectors

        m02 = x1 * y2 - x2 * y1
        m12 = x0 * y2 - x2 * y0
        m22 = x0 * y1 - x1 * y0

        a = y0 ** 2 * m02 - y1 ** 2 * m12 + y2 ** 2 * m22
        b = 2 * (x0 * y0 * m02 - x1 * y1 * m12 + x2 * y2 * m22)
        c = x0 ** 2 * m02 - x1 ** 2 * m12 + x2 ** 2 * m22

        return a, b, c

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
    def _ids_boundary(self):
        """return the edge IDs of the boundary of the hinge
        starting in the NE and moving Clockwise"""

        label_tri, label_edge = self.id_edge
        label_tri_opp, label_edge_opp = self.id_edge_opp

        SE = (label_tri, (label_edge + 1) % 3)
        NE = (label_tri, (label_edge + 2) % 3)
        NW = (label_tri_opp, (label_edge_opp + 1) % 3)
        SW = (label_tri_opp, (label_edge_opp + 2) % 3)

        return NE, SE, SW, NW

    def plot(self):
        v0, v1, v2 = self.vectors
        vertices_t1 = [sage.all.zero_vector(2), v0, v1]
        vertices_t2 = [sage.all.zero_vector(2), v1, v2]
        return sage.all.polygon2d(vertices_t1, fill=False).plot() + sage.all.polygon2d(vertices_t2, fill=False).plot()


class Triangulation(namedtuple("Triangulation", ["triangles", "gluings"])):

    @classmethod
    def _from_flatsurf(cls, X):
        DT = X.delaunay_triangulation()

        DT_polygons = [DT.polygon(i) for i in range(DT.num_polygons())]

        triangles = [Triangle(*[sage.all.vector(edge) for edge in polygon.edges()])
                     for polygon in DT_polygons]

        gluings = {edge[0]: edge[1] for edge in DT.edge_iterator(gluings=True)}

        return Triangulation(triangles, gluings)

    @classmethod
    def square_torus(cls):
        return cls._from_flatsurf(flatsurf.translation_surfaces.square_torus())

    @classmethod
    def regular_octagon(cls):
        k = QuadraticField(2)
        sqrt2 = k.gen()

        a = sage.all.vector([1, 0])
        b = sage.all.vector([1 / sqrt2, 1 / sqrt2])
        c = sage.all.vector([0, 1])
        d = sage.all.vector([-1 / sqrt2, 1 / sqrt2])

        triangles = [Triangle(a, b, -a - b),
                     Triangle(-a, -b, a + b),
                     Triangle(c, -c - b - a, a + b),
                     Triangle(-c, c + b + a, -a - b),
                     Triangle(-d, a + b + c, -c - b - a + d),
                     Triangle(d, -a - b - c, c + b + a - d)]

        gluings = {(0, 0): (1, 0), (0, 1): (1, 1), (0, 2): (2, 2), (1, 2): (3, 2),
                   (2, 0): (3, 0), (2, 1): (4, 1), (3, 1): (5, 1), (4, 0): (5, 0), (4, 2): (5, 2)}
        gluings.update({v: k for k, v in gluings.items()})

        return Triangulation(triangles, gluings)

    @classmethod
    def arnoux_yoccoz(cls, g):
        if g < 3:
            raise ValueError("g must be >= 3")
        return cls._from_flatsurf(flatsurf.translation_surfaces.arnoux_yoccoz(g))

    @classmethod
    def octagon_and_squares(cls):
        return cls._from_flatsurf(flatsurf.translation_surfaces.octagon_and_squares())

    @staticmethod
    def _triangulate_rectangle(base, height):
        triangle_lower = Triangle(sage.all.vector([0, height]),
                                  sage.all.vector([-base, -height]),
                                  sage.all.vector([base, 0]))
        triangle_upper = Triangle(sage.all.vector([0, -height]),
                                  sage.all.vector([base, height]),
                                  sage.all.vector([-base, 0]))
        return [triangle_lower, triangle_upper]

    @classmethod
    def mcmullen_l(cls, a, b):
        if not (a > 1 and b > 1):
            raise ValueError("Need to have a, b > 1")
        elif a.parent() != b.parent():
            raise ValueError("a, b need to come from same field")
        triangles = [*Triangulation._triangulate_rectangle(2, 2),
                     *Triangulation._triangulate_rectangle(b - 1, 2),
                     *Triangulation._triangulate_rectangle(2, a - 1),
                     *Triangulation._triangulate_rectangle(b - 1, 2),
                     *Triangulation._triangulate_rectangle(2, a - 1)]

        gluings = {(0, 0): (3, 0), (0, 1): (1, 1), (0, 2): (9, 2),
                   (1, 0): (6, 0), (1, 2): (4, 2), (2, 0): (7, 0),
                   (2, 1): (3, 1), (2, 2): (3, 2), (4, 0): (5, 0),
                   (4, 1): (5, 1), (5, 2): (8, 2), (6, 1): (7, 1),
                   (6, 2): (7, 2), (8, 0): (9, 0), (8, 1): (9, 1)}
        gluings.update({v: k for k, v in gluings.items()})

        return Triangulation(triangles, gluings)

    @classmethod
    def mcmullen_s(cls, a):
        triangles = [tri for dimensions in [(1, 1), (1 + a, 1), (1 + a, a), (a, a)]
                     for tri in Triangulation._triangulate_rectangle(*dimensions)]
        gluings = {(0, 0): (3, 0), (0, 1): (1, 1), (0, 2): (1, 2),
                   (1, 0): (2, 0), (3, 2): (4, 2), (2, 1): (3, 1),
                   (2, 2): (5, 2), (4, 0): (7, 0), (4, 1): (5, 1),
                   (5, 0): (6, 0), (6, 1): (7, 1), (6, 2): (7, 2)}
        gluings.update({v: k for k, v in gluings.items()})

        return Triangulation(triangles, gluings)

    @classmethod
    def ronen_l(cls, d):
        if d < 5 or d % 4 not in [0, 1]:
            raise ValueError("Must have d = 0 or 1 mod 4, d >= 5")
        c = 0 if d % 2 == 0 else -1
        length_rectangle = (d - c * c) // 4

        if not ZZ(d).is_square():
            k_rootd = QuadraticField(d)
            rootd = k_rootd.gen()
        else:
            rootd = sqrt(ZZ(d))

        side_square = (c + rootd) / 2

        dimensions = [(side_square, 1), (length_rectangle - side_square, 1), (side_square, side_square)]
        triangles = [tri for dims in dimensions for tri in Triangulation._triangulate_rectangle(*dims)]

        gluings_boundary = {(0, 2): (5, 2), (2, 2): (3, 2), (2, 0): (1, 0), (4, 0): (5, 0)}
        gluings_interior = {(0, 1): (1, 1), (2, 1): (3, 1), (4, 1): (5, 1),
                            (0, 0): (3, 0), (1, 2): (4, 2)}

        gluings = {**gluings_interior, **gluings_boundary}
        gluings.update({val: key for key, val in gluings.items()})

        return Triangulation(triangles, gluings)

    def neighbors(self, idx_tri):
        return [self.gluings[(idx_tri, k)][0] for k in range(3)]

    @property
    def edges(self):
        return [edge for edge in itertools.product(range(len(self.triangles)), range(3))
                if edge < self.gluings[edge]]

    @property
    def hinges(self):
        return [Hinge._from_id_edge(self, edge) for edge in self.edges]

    def is_delaunay(self, non_degenerate=False):
        return all(hinge.incircle_det > 0
                   if non_degenerate
                   else hinge.incircle_det >= 0
                   for hinge in self.hinges)

    def make_delaunay(self):
        while not self.is_delaunay():
            idx = randint(0, len(self.hinges) - 1)
            h = self.hinges[idx]
            if h.is_convex and h.incircle_det < 0:
                return self.flip_hinge(h.id_edge).make_delaunay()

        return self

    def change_delaunay_triangulation(self):
        if not self.is_delaunay():
            raise ValueError("Starting triangulation isn't Delaunay")
        elif self.is_delaunay(non_degenerate=True):
            return self
        while True:
            idx = randint(0, len(self.hinges) - 1)
            h = self.hinges[idx]
            if h.is_convex and h.incircle_det == 0:
                return self.flip_hinge(h.id_edge)

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
        NE, SE, SW, NW = hinge_flipped._ids_boundary

        IDs_new = {NE: SE, SE: SW, SW: NW, NW: NE}

        if edge in IDs_new:
            return IDs_new[edge]
        return edge

    def _gluings_after_flip(self, hinge_flipped):
        return {self._id_edge_after_flip(hinge_flipped, key):
                    self._id_edge_after_flip(hinge_flipped, value)
                for key, value in self.gluings.items()}

    def flip_hinge(self, id_edge):
        hinge = Hinge._from_id_edge(self, id_edge)
        if not hinge.is_convex:
            raise ValueError("Cannot flip concave hinge")
        hinge_flipped = hinge.flip()

        return Triangulation(self._triangles_after_flip(hinge_flipped),
                             self._gluings_after_flip(hinge_flipped))

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
    def idr(self):
        halfplane_to_ids_hinge = self._halfplanes_to_hinges_degen
        halfplanes = list(halfplane_to_ids_hinge.keys())

        p = halfplane.HalfPlane.intersect_halfplanes(halfplanes)

        # TODO: I'm excluding degenerate polygons as IDR due to this
        if p is None:
            return idr.IDR(p, {}, self)

        labels_segment = {idx: halfplane_to_ids_hinge[segment.halfplane]
                          for idx, segment in enumerate(p.edges)}

        return idr.IDR(p, labels_segment, self)

    def iso_delaunay_complex(self, upper_bound):
        idr_start = self.idr

        polygons_visited = {idr_start.polygon}
        segments_crossed = set()

        queue = deque([idr_start])

        while len(polygons_visited) < upper_bound:
            IDR = queue.pop()
            segments_uncrossed = [(idx, segment)
                                  for idx, segment in enumerate(IDR.polygon.edges)
                                  if segment not in segments_crossed]

            for idx_segment, segment in segments_uncrossed:
                IDR_new = IDR.cross_segment(idx_segment)
                segments_crossed |= {segment, segment.reverse()}

                if IDR_new.polygon not in polygons_visited:
                    polygons_visited.add(IDR_new.polygon)
                    queue.appendleft(IDR_new)

        return list(polygons_visited)

    @property
    def code_comb(self):
        # TODO: What's the best notion of a comb_code
        return min(comb_equiv.generate_code_marked(self, tri, edge)
                   for tri in range(len(self.triangles))
                   for edge in range(3))

    @property
    def generators_veech(self):
        return algo.generators_veech(self)

    def get_vertices_neighbor(self, vertices_tri, idx_tri, idx_edge):
        vertex_start = vertices_tri[idx_edge]
        vertex_end = vertices_tri[(idx_edge + 1) % 3]

        idx_tri_opp, idx_edge_opp = self.gluings[(idx_tri, idx_edge)]
        tri_opp = self.triangles[idx_tri_opp]
        vertices_tri_opp = {idx_edge_opp: vertex_end,
                            (idx_edge_opp + 1) % 3: vertex_start,
                            (idx_edge_opp + 2) % 3: vertex_start + tri_opp[(idx_edge_opp + 1) % 3]}
        return [vertices_tri_opp[k] for k in range(3)]

    def plot(self):
        # TODO: CLEAN + separate each part of plot into separate methods
        # TODO: How to represent edge gluings? Labels on outer edges?
        # TODO: Re-rendering when triangles intersect?

        tris_seen = {0: self.triangles[0].vertices()}
        tris_to_visit = deque([0])

        while tris_to_visit:
            idx_tri_curr = tris_to_visit.pop()
            for idx_edge in range(3):
                idx_tri_nbr, idx_edge_nbr = self.gluings[(idx_tri_curr, idx_edge)]
                if idx_tri_nbr not in tris_seen:
                    vertices_curr = tris_seen[idx_tri_curr]
                    tris_seen[idx_tri_nbr] = self.get_vertices_neighbor(vertices_curr, idx_tri_curr, idx_edge)
                    tris_to_visit.appendleft(idx_tri_nbr)

        def center(vertices):
            x = sum(vx for vx, vy in vertices) / 3
            y = sum(vy for vx, vy in vertices) / 3
            return x, y

        def midpoint(v1, v2):
            v1x, v1y = v1
            v2x, v2y = v2
            return sage.all.vector([(v1x + v2x) / 2, (v1y + v2y) / 2])

        def displacement(v1, v2):
            a, b = v1
            c, d = v2
            return 1 / 30 * sage.all.vector([b - d, c - a])

        plots_labels_tri = sum(sage.all.text(str(idx), center(vertices), fontsize=14, color='orange').plot()
                               for idx, vertices in tris_seen.items())

        plots_labels_edge = sum(sage.all.text(str(idx), midpoint(v1, v2) + displacement(v1, v2)).plot()
                                for (a, b, c) in tris_seen.values()
                                for idx, (v1, v2) in [(0, (a, b)), (1, (b, c)), (2, (c, a))])

        plots_tris = sum(sage.all.polygon2d(vertices, fill=False).plot() for vertices in tris_seen.values())

        return plots_tris + plots_labels_tri + plots_labels_edge

    @property
    def dual_graph(self):
        return Graph({idx: self.neighbors(idx) for idx in range(len(self.triangles))})


if __name__ == "__main__":
    import itertools
    import bowman.comb_equiv as ce
    import bowman.geom_equiv as ge

    # fund_dom = Triangulation.ronen_l(13).generators_veech
    # trins_ce = [idr.triangulation for idr in next(iter(fund_dom.codes_to_idrs.values()))]

    X = Triangulation.arnoux_yoccoz(3)
    idr0 = X.idr
    t0 = idr0.triangulation
    # for idx, neighbor in enumerate(idr0.neighbors):
    #     for idx1, neighbor1 in enumerate(neighbor.neighbors):
    #         t = neighbor1.triangulation
    #         print(idx, idx1, ce.gen_comb_equivs(t0, t))
    trins_ce = [t0, idr0.neighbors[1].neighbors[5].triangulation]

    for trin1, trin2 in itertools.combinations(trins_ce, r=2):
        ces = ce.gen_comb_equivs(trin1, trin2)
        print(*ces, sep='\n')
        print([ge.gen_geom_equiv(trin1, trin2, x) for x in ces])
        edges_total = list(itertools.product(range(len(trin1.triangles)), range(3)))
        codes1 = [(ce.generate_code_marked(trin1, label_tri, label_edge), (label_tri, label_edge))
                  for label_tri, label_edge in edges_total]
        codes2 = [(ce.generate_code_marked(trin2, label_tri, label_edge), (label_tri, label_edge))
                  for label_tri, label_edge in edges_total]

        _, (label_t1, label_e1) = min(codes1, key=lambda x: x[0])
        _, (label_t2, label_e2) = min(codes2, key=lambda x: x[0])

        cr1 = ce.canonical_relabel(trin1, label_t1, label_e1)
        cr2 = ce.canonical_relabel(trin2, label_t2, label_e2)

        cr1inv = {v: k for k, v in cr1.items()}
        cr2inv = {v: k for k, v in cr2.items()}

        ce_fixed = {(k, 0): cr2inv[cr1[(k, 0)]] for k in range(len(trin1.triangles))}
        for label_t_other, label_e_other in edges_total:
            cr_other = ce.canonical_relabel(trin2, label_t_other, label_e_other)
            composition = {(k, 0): cr_other[ce_fixed[(k, 0)]] for k in range(len(trin1.triangles))}
            print(composition)


