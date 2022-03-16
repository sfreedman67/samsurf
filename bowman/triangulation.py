import itertools
from collections import defaultdict, deque
from functools import lru_cache
from typing import overload

from sage.all import *
from random import randrange
# import flatsurf

from bowman import halfplane
from bowman import idr
from bowman import comb_equiv
from bowman import geom_equiv
from bowman import algo
from bowman.triangle import Triangle, is_valid_barycentric_coordinate, intersect_lines, is_point_on_line
from bowman.hinge import Hinge
from bowman.radical import Radical
from bowman.geom_equiv import gen_geom_equivs, is_cut_paste_equiv

from bowman.rational_ht_application import bicuspid_segments, segments_for_plotting
from bowman.geom_equiv import gen_geom_equivs

def return_shear_mat(dir):
    """Generate a shear that projects vector dir onto the real line, or rotates
    dir by 90 degrees if dir is verticle."""
    dir_x, dir_y = dir

    if dir_x != 0:
        return matrix([[1, 0], [(-dir_y / dir_x), 1]])
    else:
        return matrix([[0, -1], [1, 0]])


class Triangulation:
    def __init__(self, triangles=None, gluings=None):
        self.triangles = tuple(triangles) if triangles is not None else tuple()
        self.gluings = gluings if gluings is not None else {}

        self._hash = hash(self.__key())

    def __key(self):
        gluings_ordered = {(e1, e2) for e1, e2 in self.gluings.items() if e1 < e2}
        gluings_safe = tuple(sorted(gluings_ordered))
        return self.triangles, gluings_safe

    def __hash__(self):
        return self._hash

    def __eq__(self, other):
        if type(other) is type(self):
            return self.__key() == other.__key()
        else:
            return False

    def __repr__(self):
        return f"Triangulation with {len(self.triangles)} triangles"

    def copy(self):
        # create a copy which doesn't mutate the original
        return Triangulation(self.triangles, self.gluings.copy())

    @classmethod
    def _from_flatsurf(cls, trin):
        DT = trin.delaunay_triangulation()

        DT_polygons = [DT.polygon(k) for k in range(DT.num_polygons())]

        triangles = [Triangle(*[sage.all.vector(edge) for edge in x.edges()])
                     for x in DT_polygons]

        gluings = {edge[0]: edge[1] for edge in DT.edge_iterator(gluings=True)}

        return Triangulation(triangles, gluings)

    @classmethod
    def square_torus(cls):
        e0 = sage.all.vector([QQ(1), QQ(0)])
        e1 = sage.all.vector([QQ(0), QQ(1)])
        t0 = Triangle(-e0 - e1, e0, e1)
        t1 = Triangle(e0 + e1, -e0, -e1)
        gluings = {(0, 0): (1, 0), (0, 1): (1, 1), (0, 2): (1, 2)}
        gluings.update({v: k for k, v in gluings.items()})
        return Triangulation([t0, t1], gluings)

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
    def prym_eigenform_type_aplus(cls, w, h, t, e):
        """Constructs as a polygon the Type A+ Prym eigenform corresponding to
        the four nonnegative integers w, h, t, e, as detailed in the paper by
        Lanneau-Nguyen."""

        # Verify that the input satisfies the necessary conditions.
        assert w > 0 and h > 0 and t >= 0
        assert t < gcd(w, h) and gcd(gcd(gcd(w, h), t), e) == 1
        D = e ** 2 + 8 * w * h
        k = QuadraticField(D)
        sqrtD = k.gen()
        l = (e + sqrtD) / 2
        assert sign(l) > 0 and sign(w - l) > 0

        east_long = sage.all.vector(k, [w, 0])
        east_short = sage.all.vector(k, [l, 0])
        north_short = sage.all.vector(k, [0, l])
        north_east = sage.all.vector(k, [t, h])

        triangles = [Triangle(east_short, -east_short + north_short, -north_short),
                     Triangle(-east_short, east_short - north_short, north_short),
                     Triangle(north_east, -east_short, east_short - north_east),
                     Triangle(east_short, -east_short + north_east, -north_east),
                     Triangle(-east_long + east_short, east_long - east_short - north_east, north_east),
                     Triangle(east_long - east_short, -east_long + east_short + north_east, -north_east),
                     Triangle(north_east, -east_short, east_short - north_east),
                     Triangle(east_short, -east_short + north_east, -north_east),
                     Triangle(-east_long + east_short, east_long - east_short - north_east, north_east),
                     Triangle(east_long - east_short, -east_long + east_short + north_east, -north_east)]

        gluings = {(0, 0): (6, 1), (0, 1): (1, 1), (0, 2): (1, 2), (1, 0): (3, 0),
                   (2, 0): (5, 2), (2, 1): (7, 0), (2, 2): (3, 1), (3, 2): (4, 2),
                   (4, 0): (5, 0), (4, 1): (5, 1), (6, 0): (9, 2), (6, 2): (7, 1),
                   (7, 2): (8, 2), (8, 0): (9, 0), (8, 1): (9, 1)}
        gluings.update({v: k for k, v in gluings.items()})

        return Triangulation(triangles, gluings)

    @classmethod
    def prym_eigenform_type_aminus(cls, w, h, t, e):
        """Constructs as a polygon the Type A- Prym eigenform corresponding to
        the four nonnegative integers w, h, t, e, as detailed in the paper by
        Lanneau-Nguyen."""

        # Verify that the input satisfies the necessary conditions.
        assert w > 0 and h > 0 and t >= 0
        assert t < gcd(w, h) and gcd(gcd(gcd(w, h), t), e) == 1
        D = e ** 2 + 8 * w * h
        k = QuadraticField(D)
        sqrtD = k.gen()
        l = (e + sqrtD) / 2
        assert sign(l) > 0 and sign(w - l) > 0

        east_long = sage.all.vector(k, [w - l, 0])
        east_short = sage.all.vector(k, [l / 2, 0])
        north_short = sage.all.vector(k, [0, l / 2])
        north_east = sage.all.vector(k, [t, h])

        triangles = [Triangle(east_short, north_east - east_short, -north_east),
                     Triangle(-east_short, -north_east + east_short, north_east),
                     Triangle(east_short, north_east - east_short, -north_east),
                     Triangle(-east_short, -north_east + east_short, north_east),
                     Triangle(east_long, north_east - east_long, - north_east),
                     Triangle(-east_long, -north_east + east_long, north_east),
                     Triangle(east_short, north_short - east_short, -north_short),
                     Triangle(-east_short, -north_short + east_short, north_short),
                     Triangle(east_short, north_short - east_short, -north_short),
                     Triangle(-east_short, -north_short + east_short, north_short)]

        gluings = {(0, 0): (7, 0), (0, 1): (1, 1), (0, 2): (5, 2), (1, 0): (6, 0),
                   (1, 2): (2, 2), (2, 0): (9, 0), (2, 1): (3, 1), (3, 0): (8, 0),
                   (3, 2): (4, 2), (4, 0): (5, 0), (4, 1): (5, 1), (6, 1): (7, 1),
                   (6, 2): (7, 2), (8, 1): (9, 1), (8, 2): (9, 2)}
        gluings.update({v: k for k, v in gluings.items()})

        return Triangulation(triangles, gluings)

    @classmethod
    def prym_eigenform_type_b_disc_8_fake(cls):
        intersections_cylinder = Graph({0: [3, 4], 1: [3, 4, 5], 2: [4, 5]})
        perm_north = {(0, 3): (1, 3), (0, 4): (2, 4), (1, 3): (0, 3), (1, 4): (0, 4),
                      (1, 5): (2, 5), (2, 4): (1, 4), (2, 5): (1, 5)}
        perm_east = {(0, 3): (0, 4), (0, 4): (0, 3), (1, 3): (1, 4), (1, 4): (1, 5),
                     (1, 5): (1, 3), (2, 4): (2, 5), (2, 5): (2, 4)}
        # Corresponds to type B (1, 1, 0, 0)
        k = QuadraticField(8)
        a = k.gen()
        l = a / 2
        heights = [l / 2, k(1), l / 2, 1 - (l / 2), l - 1, 1 - (l / 2)]

        return Triangulation.grid_graph(intersections_cylinder, perm_north, perm_east, heights)

    @classmethod
    def prym_eigenform_type_b_disc_8_real(cls):
        k = QuadraticField(8)
        rt8 = k.gen()
        lmbd = rt8 / 2

        east_short = sage.all.vector(k, [1 - (lmbd / 2), 0])
        east_long = sage.all.vector(k, [lmbd - 1, 0])
        north_short = sage.all.vector(k, [0, lmbd / 2])
        north_long = sage.all.vector(k, [0, 1])

        t0 = Triangle(east_short, north_short, -east_short - north_short)
        t1 = Triangle(east_short + north_short, -east_short - east_long , -north_short + east_long)
        t2 = Triangle(north_short - east_long, -north_short, east_long)
        t3 = Triangle(-east_long, east_long - north_short, north_short)
        t4 = Triangle(-north_short, north_short + east_short, -east_short)
        t5 = Triangle(east_short, -east_short + north_long, -north_long)
        t6 = Triangle(north_long, -east_short, east_short - north_long)
        t7 = Triangle(north_long - east_short, -east_long - north_long, east_short + east_long)
        t8 = Triangle(east_long + north_long, -east_long - east_short, -north_long + east_short)
        t9 = Triangle(east_short + east_long, -east_long + north_short, -east_short - north_short)

        triangles = [t0, t1, t2, t3, t4, t5, t6, t7, t8, t9]

        gluings = {(0, 0): (6, 1), (0, 1): (2, 1), (0, 2): (1, 0),
                   (1, 1): (7, 2), (1, 2): (2, 0), (2, 2): (3, 0),
                   (3, 1): (9, 1), (3, 2): (4, 0), (4, 1): (9, 2),
                   (4, 2): (5, 0), (5, 1): (8, 2), (5, 2): (6, 0),
                   (6, 2): (7, 0), (7, 1): (8, 0), (8, 1): (9, 0)}
        gluings.update({v: k for k, v in gluings.items()})

        return Triangulation(triangles, gluings)

    @classmethod
    def grid_graph(cls, graph_intersections_cylinder, perm_north, perm_east, heights):
        intersections_cylinder = list(graph_intersections_cylinder.edge_iterator(labels=False))
        rectangles = [Triangulation.triangulate_rectangle(heights[idx_v], heights[idx_h])
                      for idx_h, idx_v in intersections_cylinder]
        triangles = [tri for rectangle in rectangles for tri in rectangle.triangles]
        # glue diagonal of each rectangle
        gluings = {(k, 1): (k + 1, 1) for k in range(0, len(triangles), 2)}
        for idx_intersection, intersection in enumerate(intersections_cylinder):
            # glue right edge of R_e to left edge of R_{east[e]}
            idx_intersection_east = intersections_cylinder.index(perm_east[intersection])
            gluings[(2 * idx_intersection, 0)] = (2 * idx_intersection_east + 1, 0)
            # glue top edge of R_e to bottom edge of R_{north[e]}
            idx_intersection_north = intersections_cylinder.index(perm_north[intersection])
            gluings[(2 * idx_intersection + 1, 2)] = (2 * idx_intersection_north, 2)
        gluings.update({(v, k) for k, v in gluings.items()})
        return Triangulation(triangles, gluings)

    @staticmethod
    def triangulate_rectangle(base, height):
        """ return a triangulated rectangle. Diagonal runs from lower left to upper right"""
        triangle_lower = Triangle(sage.all.vector([QQ(0), height]),
                                  sage.all.vector([-base, -height]),
                                  sage.all.vector([base, QQ(0)]))
        triangle_upper = Triangle(sage.all.vector([QQ(0), -height]),
                                  sage.all.vector([base, height]),
                                  sage.all.vector([-base, QQ(0)]))
        return Triangulation([triangle_lower, triangle_upper], {(0, 1): (1, 1), (1, 1): (0, 1)})

    @classmethod
    def mcmullen_l(cls, a, b):
        if not (a > 1 and b > 1):
            raise ValueError("Need to have a, b > 1")
        elif a.parent() != b.parent():
            raise ValueError("a, b need to come from same field")
        triangles = [*Triangulation.triangulate_rectangle(2, 2).triangles,
                     *Triangulation.triangulate_rectangle(b - 1, 2).triangles,
                     *Triangulation.triangulate_rectangle(2, a - 1).triangles,
                     *Triangulation.triangulate_rectangle(b - 1, 2).triangles,
                     *Triangulation.triangulate_rectangle(2, a - 1).triangles]

        gluings = {(0, 0): (3, 0), (0, 1): (1, 1), (0, 2): (9, 2),
                   (1, 0): (6, 0), (1, 2): (4, 2), (2, 0): (7, 0),
                   (2, 1): (3, 1), (2, 2): (3, 2), (4, 0): (5, 0),
                   (4, 1): (5, 1), (5, 2): (8, 2), (6, 1): (7, 1),
                   (6, 2): (7, 2), (8, 0): (9, 0), (8, 1): (9, 1)}
        gluings.update({v: k for k, v in gluings.items()})

        return Triangulation(triangles, gluings)

    @classmethod
    def mcmullen_s(cls, a):
        triangles = [tri for dimensions in [(QQ(1), QQ(1)), (QQ(1) + a, QQ(1)), (QQ(1) + a, a), (a, a)]
                     for tri in Triangulation.triangulate_rectangle(*dimensions).triangles]
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

        dimensions = [(side_square, k_rootd(1)), (length_rectangle - side_square, k_rootd(1)),
                      (side_square, side_square)]
        triangles = [tri for dims in dimensions for tri in Triangulation.triangulate_rectangle(*dims).triangles]

        gluings_boundary = {(0, 2): (5, 2), (2, 2): (3, 2), (2, 0): (1, 0), (4, 0): (5, 0)}
        gluings_interior = {(0, 1): (1, 1), (2, 1): (3, 1), (4, 1): (5, 1),
                            (0, 0): (3, 0), (1, 2): (4, 2)}

        gluings = {**gluings_interior, **gluings_boundary}
        gluings.update({val: key for key, val in gluings.items()})

        return Triangulation(triangles, gluings)

    def apply_matrix(self, m):
        tris_new = tuple(tri.apply_matrix(m) for tri in self.triangles)
        return Triangulation(tris_new, self.gluings)

    def apply_gt_flow(self, t):
        """Apply g_t flow to the triangulation for a time t.
        (Positive values of t coorespondes to contraction in the x-direction)"""
        g_t = matrix([[ZZ(2) ** (-t), 0], [0, ZZ(2) ** t]])
        return self.apply_matrix(g_t)

    @property
    def self_geom_equivs(self):
        return gen_geom_equivs(self, self)

    def geom_equiv_relabelling(self, equiv_trin, tri_idx, edge_idx=None):
        """
        Given a cut-and-paste equivalent triangulation to self, tells you what
        triangle or edge in self corresponds to the input indices on the other
        triangulation
        INPUT:
        * equiv_trin: a triangulation, which is cut-and-paste equivalent to self
        (i.e. geometric equivalent with the identity as the equivalence matrix)
        * tri_idx: int, index of triangle on equiv_trin
        * edge_idx: if None, returns triangle on self corresponding to tri_idx
                    if int in [0, 1, 2], returns edge on self corresponding to
                    (tri_idx, edge_idx) on equiv_trin
        """
        equivs = gen_geom_equivs(equiv_trin, self)
        # geometric equivalences from equiv_trin to self
        equiv = [e for e in equivs if e[0] == matrix([[1, 0], [0, 1]])][0]
        # a cut and paste equivalence between the two
        edge_rel_dict = equiv[1]  # the edge relabelings stored in equiv[1]
        tri_rel_dict = {e1[0]: e2[0] for e1, e2 in edge_rel_dict.items()}
        if edge_idx is None:
            return tri_rel_dict[tri_idx]
        elif edge_idx in [0, 1, 2]:
            return edge_rel_dict[(tri_idx, edge_idx)]
        else:
            error_msg = (f"edge_idx must be None for triangle relabelling or" +
                         f" in [0, 1, 2] for edge relabelling, not {edge_idx}")
            raise ValueError(error_msg)

    def geom_equiv_relabel_point(self, equiv_trin, tri_idx, pt_coords):
        """
        Given a cut-and-paste equivalent triangulation to self, finds the new
        location of a given point (triangle id + bary coordinates inside that
        triangle)
        OUTPUT:
        tuple (new_tri_idx, new_pt_coords)
        """
        new_tri_idx = self.geom_equiv_relabelling(equiv_trin, tri_idx)
        edge_relabelling = {
            edge: self.geom_equiv_relabelling(equiv_trin, tri_idx, edge)[1]
            for edge in [0, 1, 2]}  # how edges of the triangle are permuted
        new_pt_coords = [None] * 3  # initialize coordinates
        for edge_id, relabelling in edge_relabelling.items():
            new_pt_coords[relabelling - 1 % 3] = pt_coords[edge_id - 1 % 3]
            # barycentric coordinates (x:y:z) correspond to heights from edges
            # (1, 2, 0), so we offset by 1 mod 3
            # e.g. if relabelling of 2 is 0, then midpoint of 2 gets sent to
            # midpoint of 0
        return (new_tri_idx, tuple(new_pt_coords))

    def check_horiz(self):
        """Check if every triangle has a horizontal edge.  Returns True if so."""
        tris = self.triangles
        v = vector([1, 0])
        for i in tris:
            horiz = False
            for j in range(0, 3):
                if (abs(v.dot_product(i[j])) == i[j].norm() * v.norm()):
                    horiz = True
            if (horiz):
                continue
            else:
                return False
        return True

    def get_horiz_cyls(self):
        """Take triangulation with every triangle having a horizontal edge,
        and return the horizontal cylinder refinement."""
        tris = list(self.triangles)
        gluings = self.gluings
        v = vector([1, 0])
        cylinders = []
        while True:
            if (not tris):
                return cylinders
            # make cylinder
            cylinder = []
            curr_tri = tris.pop(0)
            while True:
                curr_tri_indx = self.triangles.index(curr_tri)
                cylinder.append(curr_tri_indx)
                check = False
                for i in range(0, 3):
                    if (abs(v.dot_product(curr_tri[i])) != curr_tri[i].norm()):
                        glued_tri = self.triangles[gluings[(curr_tri_indx, i)][0]]
                        glued_tri_indx = self.triangles.index(glued_tri)
                        if ((glued_tri_indx not in cylinder) and (glued_tri in tris)):
                            tris.remove(glued_tri)
                            check = True
                            curr_tri = glued_tri
                            break
                if (not check):
                    break
            cylinders.append(cylinder)
        return cylinders

    def make_directional_triangulation(self, direction):
        """This function takes a triangulation of a translation surface along with a
        cylinder direction, and returns a triangulation where all triangles have
        an edge parallel to the specified cylinder direction.  It also returns the cylinder
        refinement for each cylinder of the surface, in the specified direction.

        INPUT.      direction := a 2D vector pointing in the cylinder direction of interest.
        OUTPUT.     new_triangulation := triangulation with every triangle having edge in
                        specified direction
                    cylinders := a list, each of whose elements are lists of indices corresponding
                        to the triangles in the refinement of a cylinder.
        """
        mat = return_shear_mat(direction)
        matinv = mat.inverse()
        new_triangulation = self.apply_matrix(mat)
        counter = 0

        while True:
            if(new_triangulation.is_delaunay):
                if(new_triangulation.check_horiz()):
                    # print("Completed triangulation.")
                    break
            # else apply g_t flow until no-longer delaunay, and retriangulate
            while True:
                if (new_triangulation.is_delaunay):
                    counter += 1
                    new_triangulation = new_triangulation.apply_gt_flow(2)
                else:
                    new_triangulation = new_triangulation.make_delaunay()
                    break
            if (counter >= 25):
                print("Exited loop after applying g_t flow for 25 iterations")
                break
        # get cylinders, then g_t flow in inverse direction and rotate back
        cylinders = new_triangulation.get_horiz_cyls()
        new_triangulation = new_triangulation.apply_gt_flow(-(counter * 2))
        new_triangulation = new_triangulation.apply_matrix(matinv)
        return new_triangulation, cylinders

    def compute_constraints_transformed(self, veech_elem):
        """Returns a list of line segments (pairs of barycentric coordinates)
        that are defined by the constraints given by the Rational Height Lemma.
        VEECH_ELEM transforms these line segments before their return."""
        constraints_dict = segments_for_plotting(bicuspid_segments(self))
        marked_tris = self.plot_constraints(constraints_dict).main_constraint_plotter(veech_elem).triangles
        # only add the transformed constraints, not including the originals (distinguished by color)
        constraints_transformed_list = [[(line[0],line[1]) for line in marked_tris[i].lines_marked if line[2] == (0,0.9,0.1)] for i in range(len(self.triangles))]
        return constraints_transformed_list

    def compute_candidate_periodic_points(self, tri_id, veech_elem):
        # Input := tri_id for triangle in self, and veech_elem to apply to the constraints
        #       in the triangle, and then compute the intersections.  Returns intersection points.

        # 1. Identify those generators for which the transformed constraints
        #    can be computed.
        # good_gens = veech_gens_list
        # good_gens = good_gens + [gen**(-1) for gen in good_gens]
        #num_good_gens = len(good_gens)
        #print(f"Identified {num_good_gens} good generators.")

        # Obtain initial constraints
        gen = matrix([[1,0],[0,1]])
        initial_constraints = self.compute_constraints_transformed(gen)[tri_id]
        lines = {geo_elem for geo_elem in initial_constraints if geo_elem[0] != geo_elem[1]}
        points = {geo_elem[0] for geo_elem in initial_constraints if geo_elem[0] == geo_elem[1]}

        while lines:
            #print(f"Number of lines to eliminate: {len(lines)}.")
            # print("Lines:")
            # for line in lines:
            #     print(line)

            # OLD CODE that randomly applied elements from generators
            # gen = gen * good_gens[randrange(num_good_gens)]
            gen = gen * veech_elem

            #print(f"Applying {gen}...")
            new_constraints = self.compute_constraints_transformed(gen)[tri_id]
            new_lines = {geo_elem for geo_elem in new_constraints if geo_elem[0] != geo_elem[1]}
            new_points = {geo_elem[0] for geo_elem in new_constraints if geo_elem[0] == geo_elem[1]}

            # a. Intersect the current points with...
            # ...the new points.
            next_points = list(points & new_points)
            # ...the new lines.
            for line in new_lines:
                for point in points:
                    if is_point_on_line(point, line):
                        next_points.append(point)

            # b. Intersect the current lines with the new points.
            for line in lines:
                for point in new_points:
                    if is_point_on_line(point, line):
                        next_points.append(point)

            # c. Intersect the current lines with the new lines.
            next_lines = []
            for line in lines:
                for new_line in new_lines:
                    is_intersect, local_obj = intersect_lines(*line, *new_line)
                    if is_intersect:
                        if len(local_obj) == 2:
                            next_lines.append(local_obj)
                        else:
                            next_points.append(local_obj)

            points = set(next_points)
            lines = set(next_lines)
        return points

    def _return_triangle_coords(self, id_tri):
        """Helper funciton returning triangle orientation type for rectilinear triangle
        located at id_tri."""
        tri = self.triangles[id_tri]
        for k, edge in enumerate(tri):
            x, y = edge
            if (y == 0):
                vec_horiz = edge
                idx_horiz = k
            elif (x == 0):
                vec_vert = edge
                idx_vert = k
            else:
                vec_hyp = edge
                idx_hyp = k

        # check orientation (with right angle at the origin, type := quadrant containing the hypotenuse)
        if (vec_vert[1] < 0 and vec_horiz[0] > 0):
            # type 1
            coords_horiz = vector([QQ(0), QQ(0)])
            coords_vert = -vec_vert
            coords_hyp = vec_horiz
        elif (vec_vert[1] > 0 and vec_horiz[0] > 0):
            # type 2
            coords_horiz = vector([QQ(0), QQ(0)])
            coords_vert = vec_horiz
            coords_hyp = vec_horiz + vec_vert
        elif (vec_vert[1] > 0 and vec_horiz[0] < 0):
            # type 3
            coords_horiz = -vec_horiz + vec_vert
            coords_vert = -vec_horiz
            coords_hyp = vec_vert
        else:
            # type 4
            coords_horiz = vec_hyp
            coords_vert = -vec_vert
            coords_hyp = vector([QQ(0), QQ(0)])

        idx_to_coords = {idx_vert: coords_vert, idx_horiz: coords_horiz, idx_hyp: coords_hyp}
        return [idx_to_coords[k] for k in range(3)]

    def _return_line_params(self, p1, p2):
        """Take two cartesian coordinates p1 and p2, and return the two parameters
        which determine the line cutting through both."""
        r = p1[1] - p2[1]
        ru = p1[0] - p2[0]
        m = 0
        if (ru == 0):
            return (oo, p1[0])
        else:
            m = r / ru  # need QQ?
        b = p1[1] - (m * p1[0])
        return (m, b)

    def _get_intersection_pt(self, line1, line2):
        """Takes two line parameters and returns intersection pt in
        cartesian coordinates.  line1/2 = (m, b) for slope m, intercept b."""
        if (line1[0] == oo):
            return (line1[1], (line1[1] * line2[0]) + line2[1])
        if (line2[0] == oo):
            return (line2[1], (line2[1] * line1[0]) + line1[1])
        if (line1[0] == line2[0]):
            raise ValueError("Lines are parallel.")
        x = -(line1[1] - line2[1]) / (line1[0] - line2[0])
        y = line1[0] * x + line1[1]
        return (x, y)

    def return_intersections(self, triangle_id, constraints_list):
        """Take a triangle and a list of constraints of form [a, b, c]
        where a, b, c are the parameters that define a line (a for x, b for y).
        Return barycentric coordinates that correspond to where on the triangle
        the lines intersect, with respect to the determined coordinate system.

        triangle_id         := the index of the triangle in self.triangles.
        constraints list    := list of constraints, each constrain is a list of three parameters
        RETURNS:            collection of points, corresponding to the points at which 
                            each constraint passes through the triangle edges.
        """
        # unpack coordinates of triangle with appropriate origin location.
        r = self._return_triangle_coords(triangle_id)
        # get triangle lines
        tri_lines = []
        for i in range(3):
            j = (i + 1) % 3
            p1, p2 = r[i], r[j]
            line = self._return_line_params(p1, p2)
            tri_lines.append(line)
        # check where constraint lines intersect triangle lines.
        # make sure not to check parallel lines.
        line_intersections = []
        for i in range(len(constraints_list)):
            # make sure not dividing by 0
            # cons_line = (0, 0)
            if (constraints_list[i][1] == 0):
                cons_line = (oo, -constraints_list[i][2] / constraints_list[i][0])
            else:
                line_slope = -constraints_list[i][0] / constraints_list[i][1]
                line_intercept = -constraints_list[i][2] / constraints_list[i][1]
                cons_line = (line_slope, line_intercept)
            intersections = []
            for j in range(3):
                tri_line = tri_lines[j]
                if (tri_line[0] == cons_line[0]):  # parallel
                    continue
                intersection_pt = self._get_intersection_pt(cons_line, tri_line)
                intersections.append(intersection_pt)
            line_intersections.append(intersections)
        # Keep the points that lie on a triangle edge.
        x_min = min([r[0][0], r[1][0], r[2][0]])
        x_max = max([r[0][0], r[1][0], r[2][0]])
        y_min = min([r[0][1], r[1][1], r[2][1]])
        y_max = max([r[0][1], r[1][1], r[2][1]])
        kept_coords = []
        for i in range(len(line_intersections)):
            good_coords = []
            for j in range(len(line_intersections[i])):
                pt = line_intersections[i][j]
                if (pt[0] < x_min or pt[0] > x_max or pt[1] < y_min or pt[1] > y_max):
                    continue
                else:
                    pt = vector(pt)
                    good_coords.append(pt)
            kept_coords.append(good_coords)

        # convert intersection pts to barycentric coordinates
        return [[Triangulation.cart_to_bary(p, r) for p in coords] for coords in kept_coords]

    @staticmethod
    def cart_to_bary(cart_pt, r):
        """
            cart_py := cartesian coordinate
            r := list of three triangle vertices
        """
        (x0, y0), (x1, y1), (x2, y2) = r
        x, y = cart_pt

        T = matrix([[x0 - x2, x1 - x2, x2],
                    [y0 - y2, y1 - y2, y2],
                    [0, 0, 1]])

        lambda1, lambda2, _ = T.inverse() * vector([x, y, 1])
        return (lambda1, lambda2, 1 - lambda1 - lambda2)

    def plot_constraints(self, dict):
        """Takes as input a dictionary relating each triangle to a list of constraints.
        This function then plots these constraints on the triangulation."""
        tris = self.triangles
        num_tris = len(tris)
        new_triangulation = Triangulation(self.triangles, self.gluings)
        for i in range(num_tris):
            coords = self.return_intersections(i, dict[i])  # list of pairs of barycentric coords
            for j in range(len(coords)):
                if (len(coords) == 0 or len(coords[j]) < 2):
                    continue
                else:
                    new_triangulation = new_triangulation.mark_line(i, coords[j][0], coords[j][1], (1, 0, 0))
        return new_triangulation

    def track_marked_point(self, coord, triangle_id, veech_elem):
        """
        coord := tuple of three nonnegative real numbers summing to 1
        triangle_id := triangle index of triangulation
        veech_elem := element of veech group for t the surface.
        The function marks the triangulation with the point,
        and applies the veech element.  It then returns the triangle id and new barycentric
        coordinate of the marked point, of form (tri_id, coord).
        """
        new_triangulation = Triangulation(self.triangles, self.gluings)
        new_triangulation = new_triangulation.mark_point(triangle_id, coord, (1, 0, 0))
        new_triangulation = new_triangulation.apply_matrix(veech_elem)
        new_triangulation = new_triangulation.make_delaunay(self)

        # now find the marked point and return the coordinate and triangle
        tris = new_triangulation.triangles
        for i, tri in enumerate(tris):
            pts_marked = tri.points_marked  # list containing pairs ((a, b, c), (r, g, b)) for bary and rgb color
            if pts_marked:
                real_tri_indx, real_pt = self.geom_equiv_relabel_point(new_triangulation, i, pts_marked[0][0])
                return real_tri_indx, real_pt
        return triangle_id, coord

    def plot_transformed_constraints(self, pts_info):
        """
        self     := segments triangulation
        pts_info := pts_info is a list.  Each element is a list containing a collection of triples
                of form (bary_coord, tri_indx, vector) which posses information to be passed
                into mark_flow to plot a line, starting at bary_coord, in tri_indx, in direction
                of vector.
        """
        tris = self.triangles # use self.triangles as this was order pts_info given in
        #new_triangulation = self
        new_triangulation = Triangulation(self.triangles, self.gluings)
        for tri_dat in pts_info:
            for base_pt, indx, vec in tri_dat:
                new_triangulation = new_triangulation.mark_flow(indx, base_pt, vec, 1, (0, 0.9, 0.1))
        return new_triangulation

    def main_constraint_plotter(self, veech_elem):
        """
        self   := a triangulation (with marked segments)
        veech_elem              := global veech element for triangulation
        This function returns a new triangulation with transformed segments
        under the veech element, ALONG WITH the original segments.
        """

        def midpoint_barys(p0, p1):
            return tuple((a + b) / QQ(2) for a, b in zip(p0, p1))

        def subdivide_line_marked(line_marked):
            start_orig, end_orig, color = line_marked
            midpoint_line_marked = midpoint_barys(start_orig, end_orig)
            s0 = (midpoint_line_marked, start_orig, color)
            s1 = (midpoint_line_marked, end_orig, color)
            return [s0, s1]

        new_pts_info = []
        for i, tri in enumerate(self.triangles):
            mp_info = []
            lines_subdivided = [segment for line_marked in tri.lines_marked
                                for segment in subdivide_line_marked(line_marked)]
            for base_coord, dir_coord, color in lines_subdivided:
                    if base_coord != dir_coord:
                        vector_orig = Triangulation.bary_coords_vec(base_coord, dir_coord, tri)
                        vector_transformed = veech_elem * vector_orig
                        real_tri_indx, new_coord = self.track_marked_point(base_coord, i, veech_elem)
                        mp_info.append((new_coord, real_tri_indx, vector_transformed))

            new_pts_info.append(mp_info)
        # return triangulation with BOTH original line segments AND transformed line segements
        return self.plot_transformed_constraints(new_pts_info)

    @staticmethod
    def bary_coords_vec(coord1, coord2, triangle):
        # given the two barycentric coordinates, outputs vector from one to the other
        a1, b1, c1 = coord1
        a2, b2, c2 = coord2
        v0, v1, v2 = triangle[0], triangle[1], triangle[2]
        vec1 = (b1 * v0) - (c1 * v2)
        vec2 = (b2 * v0) - (c2 * v2)
        return vec2 - vec1

    def mark_point(self, triangle_id, coords, rgbcolor):
        """ in color RGBCOLOR the point determined by barycentric coordinates COORDS on the triangle TRIANGLE_ID.

        triangle_id := the ID, i.e., index in self.triangles, of the triangle to be marked.
        coords      := a tuple of three nonnegative real numbers that sum to 1 representing the barycentric
                       coordinates of a point in the triangle determined by TRIANGLE_ID.
        rgbcolor    := a tuple of three real numbers between 0 and 1, where the componenets are the RGB values of a
                       color.
        """
        if not is_valid_barycentric_coordinate(*coords):
            raise ValueError(f"Invalid barycentric coordinates {coords}.")

        tris_new = self.triangles[0:]

        def __replace__(array, index, new_element):
            return array[:index] + (new_element,) + array[index + 1:]

        if coords[0] > 0 and coords[1] > 0 and coords[2] > 0:
            # If the point marked lies in the interior of the triangle...
            tris_new = __replace__(tris_new, triangle_id, tris_new[triangle_id].mark_point(coords, rgbcolor))
        elif coords[0] + coords[1] > 0 and coords[1] + coords[2] > 0 and coords[2] + coords[0] > 0:
            # If the point marked lies on the interior of an edge of the triangle...
            tris_new = __replace__(tris_new, triangle_id, tris_new[triangle_id].mark_point(coords, rgbcolor))
            edge_id = (coords.index(0) + 1) % 3
            opp_triangle_id, opp_edge_id = self.gluings[(triangle_id, edge_id)]
            opp_coords_indexed = sorted(
                [(opp_edge_id, coords[(edge_id + 1) % 3]), ((opp_edge_id + 1) % 3, coords[edge_id]),
                 ((opp_edge_id + 2) % 3, coords[(edge_id + 2) % 3])])
            opp_coords = tuple(opp_coord for _, opp_coord in opp_coords_indexed)
            tris_new = __replace__(tris_new, opp_triangle_id,
                                   tris_new[opp_triangle_id].mark_point(opp_coords, rgbcolor))

        return Triangulation(tris_new, self.gluings)

    def mark_line(self, triangle_id, start_coords, end_coords, rgbcolor):
        """Mark in color RGBCOLOR the line determined by barycentric coordinates START_COORDS and END_COORDS on the triangle TRIANGLE_ID.

        triangle_id  := the ID, i.e., index in self.triangles, of the triangle
                        to be marked.
        start_coords := a tuple of three nonnegative real numbers that sum to 1
                        representing the barycentric coordinates of a point in
                        the triangle determined by TRIANGLE_ID.
        rgbcolor     := a tuple of three real numbers between 0 and 1, where
                        the componenets are the RGB values of a color.
        """
        if not is_valid_barycentric_coordinate(*start_coords):
            raise ValueError(f"Invalid barycentric coordinates {start_coords}.")

        if not is_valid_barycentric_coordinate(*end_coords):
            raise ValueError(f"Invalid barycentric coordinates {end_coords}.")

        tri_new = self.triangles[triangle_id].mark_line(start_coords, end_coords, rgbcolor)
        tris_new = self.triangles[:triangle_id] + (tri_new,) + self.triangles[triangle_id + 1:]

        return Triangulation(tris_new, self.gluings)

    def __step_flow_helper__(self, start_tri_id, start_coords, direction):
        """Helper function for Triangulation.step_flow. Note that this method
        functions as expected unless DIRECTION is in the direction of a vertex
        of the starting triangle, so verify using Triangle.is_toward_conepoint
        that this is not the case before calling."""
        start_tri = self.triangles[start_tri_id]
        start_pos = start_coords[1] * start_tri[0] - start_coords[2] * start_tri[2]

        for i in range(3):
            if start_coords[i] == 1:
                assert (start_tri.is_interior(i, direction))  # Error means pointing away from triangle.

        # Step 1: Identify the outgoing edge.
        p = (start_pos, start_pos - start_tri[0], start_pos + start_tri[2])

        out_edge_is_assigned = False
        for i in range(3):
            change_of_basis = sage.all.column_matrix((-p[i], -p[(i + 1) % 3]))
            if not change_of_basis.is_singular():
                sector_coords = change_of_basis ** (-1) * direction
                if sector_coords[0].sign() > 0 and sector_coords[1].sign() > 0:
                    out_edge = i
                    out_edge_is_assigned = True

        if not out_edge_is_assigned:
            for i in range(3):
                change_of_basis = sage.all.column_matrix((-p[i], -p[(i + 1) % 3]))
                if change_of_basis.is_singular():
                    return start_coords[i], i, 0

        assert (out_edge_is_assigned)

        # Step 2: Determine the vector coordinates of the next point.
        linear_system = sage.all.column_matrix((start_tri[out_edge], direction))
        # print("martrix: ", linear_system)
        line = start_pos
        for i in range(out_edge):
            line = line - start_tri[i]
        s = (linear_system ** (-1) * line)[0]
        t = -(linear_system ** (-1) * line)[1]
        return s, out_edge, t

    def step_flow(self, start_tri_id, start_coords, direction):
        """Flows a given point in a given direction.

        start_tri_id := the triangle to which the point belongs, provided as an integer.
        start_coords := the vector from the zeroth vertex to the point.
        direction    := the direction of the flow."""

        if self.triangles[start_tri_id].is_toward_conepoint(start_coords, direction):
            return start_tri_id, start_coords, direction

        s, out_edge, _ = self.__step_flow_helper__(start_tri_id, start_coords, direction)

        # Step 3: Translate these coordinates to those of the next triangle over.
        end_tri_id, in_edge = self.gluings[(start_tri_id, out_edge)]
        end_tri = self.triangles[end_tri_id]
        end_coords_indexed = sorted([(in_edge, s),
                                     ((in_edge + 1) % 3, 1 - s),
                                     ((in_edge + 2) % 3, 0)])
        end_coords = tuple(coord for _, coord in end_coords_indexed)

        return end_tri_id, end_coords, direction

    def mark_orbit(self, start_tri_id, start_coords, direction, rgbcolor):
        """Marks the trajectory of a given point under the straight line flow
           in a given direction.

        start_tri_id := the triangle to which the point belongs, provided as an integer.
        start_coords := the vector from the zeroth vertex to the point.
        direction    := the direction of the flow."""

        points_seen = []
        tris_new = self.triangles

        while (start_tri_id, start_coords) not in points_seen:
            points_seen.append((start_tri_id, start_coords))

            # Extend the straight line from the previous point.
            assert (not self.triangles[start_tri_id].is_toward_conepoint(start_coords, direction))
            s, out_edge, t = self.__step_flow_helper__(start_tri_id, start_coords, direction)
            end_coords_indexed = sorted([((out_edge + 1) % 3, s),
                                         (out_edge, 1 - s),
                                         ((out_edge + 2) % 3, 0)])
            end_coords = tuple(coord for _, coord in end_coords_indexed)

            tris_new = tris_new[0:start_tri_id] + (
                tris_new[start_tri_id].mark_line(start_coords, end_coords, rgbcolor),) + tris_new[start_tri_id + 1:]

            # Prepare for the next depth of recursion.
            start_tri_id, in_edge = self.gluings[(start_tri_id, out_edge)]
            start_coords_indexed = sorted([(in_edge, s),
                                           ((in_edge + 1) % 3, 1 - s),
                                           ((in_edge + 2) % 3, 0)])
            start_coords = tuple(coord for _, coord in start_coords_indexed)

        return Triangulation(tris_new, self.gluings)

    def mark_flow(self, start_tri_id, start_coords, velocity, time, rgbcolor):
        """Marks the trajectory of a given point under the straight line flow
           in a given direction.

        start_tri_id := the triangle to which the point belongs, provided as an integer.
        start_coords := the vector from the zeroth vertex to the point.
        velocity     := the velocity of the flow.
        time         := the amount of time to flow for."""

        points_seen = []
        tris_new = self.triangles
        time_traveled = 0

        while (start_tri_id, start_coords) not in points_seen:
            start_tri = tris_new[start_tri_id]
            points_seen.append((start_tri_id, start_coords))
            vertex_id = sum(tuple(m + 1 for m in range(3) if start_coords[m] == 1) + (-1,))

            if tris_new[start_tri_id].is_toward_conepoint(start_coords, velocity):
                t = time - time_traveled
            elif vertex_id > -1 and not start_tri.is_interior(vertex_id, velocity):
                # TODO: Handle bad input.
                print("here")
                break
            else:
                # Otherwise, we are safe to call __step_flow_helper__.
                s, out_edge, t = self.__step_flow_helper__(start_tri_id, start_coords, velocity)
                end_coords_indexed = sorted([((out_edge + 1) % 3, s),
                                             (out_edge, 1 - s),
                                             ((out_edge + 2) % 3, 0)])
                end_coords = tuple(coord for _, coord in end_coords_indexed)

            # Figure out whether we have run the length of our trajectory.
            if time_traveled + t < time:
                tris_new = tris_new[0:start_tri_id] + (
                    start_tri.mark_line(start_coords, end_coords, rgbcolor),) + tris_new[start_tri_id + 1:]

                # Prepare for the next depth of recursion.
                start_tri_id, in_edge = self.gluings[(start_tri_id, out_edge)]
                start_coords_indexed = sorted([(in_edge, s),
                                               ((in_edge + 1) % 3, 1 - s),
                                               ((in_edge + 2) % 3, 0)])
                start_coords = tuple(coord for _, coord in start_coords_indexed)
                time_traveled = time_traveled + t
            else:
                remainder = time - time_traveled
                change_of_basis = sage.all.matrix([
                    [start_tri[0][0], -start_tri[2][0]],
                    [start_tri[0][1], -start_tri[2][1]]
                ])
                end_pos = start_coords[1] * start_tri[0] - start_coords[2] * start_tri[2] + remainder * velocity
                end_coords_partial = change_of_basis ** (-1) * end_pos
                end_coords = (
                    1 - end_coords_partial[0] - end_coords_partial[1], end_coords_partial[0], end_coords_partial[1])
                tris_new = tris_new[0:start_tri_id] + (
                    tris_new[start_tri_id].mark_line(start_coords, end_coords, rgbcolor),) + tris_new[start_tri_id + 1:]
                break

        return Triangulation(tris_new, self.gluings)


    def is_on_same_geodesic(self, start_tri_id, start_coords, end_tri_id, end_coords, direction):
        for i in range(100):
            if start_tri_id == end_tri_id:
                tri = self.triangles[start_tri_id]
                start_pos = start_coords[1] * tri[0] - start_coords[2] * tri[2]
                end_pos = end_coords[1] * tri[0] - end_coords[2] * tri[2]
                nonzero_component = 0
                if direction[nonzero_component].is_zero():
                    nonzero_component = 1
                displacement = start_pos - end_pos
                ratio = displacement[nonzero_component] / direction[nonzero_component]
                if displacement == ratio * direction:
                    return True

            start_tri_id, start_coords, direction = self.step_flow(start_tri_id, start_coords, direction)
        return False

    @property
    def edges(self):
        return [edge for edge in itertools.product(range(len(self.triangles)), range(3))
                if edge < self.gluings[edge]]

    @property
    def hinges(self):
        return [Hinge.from_id_edge(self, edge) for edge in self.edges]

    @property
    def is_delaunay(self):
        return all(hinge.incircle_det >= 0 for hinge in self.hinges)

    @property
    def is_delaunay_strict(self):
        return all(hinge.incircle_det > 0 for hinge in self.hinges)

    def make_delaunay(self, equiv_trin=None):
        """
        Makes a new Delaunay triangulation from self.
        Randomly flips hinges until it is delaunay.
        Input:
        * equiv_trin: Triangulation
            Tries to return a delaunay triangulation which is
            cut&paste equivalent to the triangulation equiv_trin
        """

        # first, find a candidate delaunay triangulation from self
        while not self.is_delaunay:
            idx = randint(0, len(self.hinges) - 1)   ###### WHERE PROBLEM IS #####
            h = self.hinges[idx]
            if h.is_convex and h.incircle_det < 0:
                self = self.flip_hinge(h.id_edge)

        candidate = self
        # now find the presentation of candidate which is cut&paste equivalent
        if equiv_trin is None:
            return candidate
        else:
            concyclic_hinges = [h.id_edge for h in candidate.hinges if
                                h.is_convex and h.incircle_det == 0]
            for possible_trin in candidate.flips_generator(concyclic_hinges):
                # checks all possible hinge flips for concyclic hinges
                if is_cut_paste_equiv(equiv_trin, possible_trin):
                    return possible_trin
            # if for loop exited, no possible_trin matches
            raise ValueError("equiv_trin has no Delaunay equivalent to self")

    def make_nontrivial(self):
        shear = sage.all.matrix([[QQ(1), QQ(0.1)], [QQ(0), QQ(1)]])
        t = self.apply_matrix(shear)
        t = t.make_delaunay()
        if not t.is_delaunay_strict:
            raise ValueError("You picked an unlucky shear!")
        return t.apply_matrix(shear.inverse())

    @property
    def halfplanes(self):
        return list(filter(lambda x: x is not None,
                           (hinge.halfplane for hinge in self.hinges)))

    @property
    def _halfplanes_to_hinges_degenerate(self):
        halfplanes_labelled = [(hinge.halfplane, hinge.id_edge)
                               for hinge in self.hinges]

        dd = defaultdict(list)
        for x, id_hinge in halfplanes_labelled:
            dd[x].append(id_hinge)

        # remove degenerate halfplanes that are None
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

    @staticmethod
    def _id_edge_after_flip(hinge_flipped, edge):
        NE, SE, SW, NW = hinge_flipped.ids_boundary

        IDs_new = {NE: SE, SE: SW, SW: NW, NW: NE}

        if edge in IDs_new:
            return IDs_new[edge]
        return edge

    def _gluings_after_flip(self, hinge_flipped):
        return {self._id_edge_after_flip(hinge_flipped, key): self._id_edge_after_flip(hinge_flipped, value)
                for key, value in self.gluings.items()}

    def flip_hinge(self, id_edge):
        hinge = Hinge.from_id_edge(self, id_edge)
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

    def flips_generator(self, ids_list):
        """
        Given a list of hinge ids, provides a generator
        Outputs of the generator are triangulations, with a subset of
        the hinge flips in ids_list applied to it
        """
        if len(ids_list) == 0:
            yield self
            return
        else:
            yield from self.flips_generator(ids_list[:-1])
            yield from [trin.flip_hinge(ids_list[-1]) for trin
                        in self.flips_generator(ids_list[:-1])]


    def plot_halfplanes(self, count=None):
        figure = sum(itertools.islice((x.plot()
                                       for x in self.halfplanes), count))
        if count is not None:
            plt_final = figure[-1]
            opt = plt_final.options()
            opt["linestyle"] = "--"
            plt_final.set_options(opt)

        return figure

    @property
    @lru_cache(None)
    def idr(self):
        halfplane_to_ids_hinge = self._halfplanes_to_hinges_degenerate
        halfplanes = list(halfplane_to_ids_hinge.keys())
        p = halfplane.HalfPlane.intersect_halfplanes(halfplanes)
        labels_segment = {idx: halfplane_to_ids_hinge[segment.halfplane]
                          for idx, segment in enumerate(p.edges)}

        return idr.IDR(p, labels_segment, self)

    def iso_delaunay_complex(self, upper_bound):
        idr_start = self.idr

        idrs_visited = {idr_start}
        segments_crossed = set()

        queue = deque([idr_start])

        while len(idrs_visited) < upper_bound:
            IDR = queue.pop()
            segments_uncrossed = [(idx, segment)
                                  for idx, segment in enumerate(IDR.polygon.edges)
                                  if segment not in segments_crossed]

            for idx_segment, segment in segments_uncrossed:
                idr_new = IDR.get_idr_neighboring(idx_segment)
                segments_crossed |= {segment, segment.reverse()}

                if idr_new.polygon not in idrs_visited:
                    idrs_visited.add(idr_new)
                    queue.appendleft(idr_new)

        return idrs_visited

    @property
    def code_comb(self):
        return next(iter(self.codes_comb))[0]

    # TODO: model as dictionary (tri, edge) --> code_comb?
    @property
    @lru_cache(None)
    def codes_comb(self):
        codes = {comb_equiv.generate_code_marked(self, tri, edge)
                 for tri in range(len(self.triangles))
                 for edge in range(3)}
        h_min = min(code[0] for code in codes)
        # minimum hash same for two comb  equiv triangulations
        return {code for code in codes if code[0] == h_min}

    @property
    def code_geom(self):
        return next(iter(self.codes_geom))[0]

    @property
    @lru_cache(None)
    def codes_geom(self):
        codes = {geom_equiv.generate_code_marked(self, edge[0], edge[1])
                 for _, edge in self.codes_comb}
        h_min = min(code[0] for code in codes)
        return {code for code in codes if code[0] == h_min}

    @property
    @lru_cache(None)
    def code(self):
        return self.code_comb, self.code_geom

    @property
    def generators_veech(self):
        return algo.generators_veech(self)

    @property
    def cylinder_directions(self):
        dir_list = list()

        for idr in self.generators_veech.idrs:
            for vertex in idr.polygon.vertices:
                if vertex == oo and (1, 0) not in [dir for dir, _ in dir_list]:
                    dir_list.append(((1, 0), idr.triangulation))
                elif vertex != oo and vertex.v2 == 0 and (vertex.u, 1) not in [dir for dir, _ in dir_list]:
                    dir_list.append(((vertex.u, 1), idr.triangulation))

        return dir_list

    def get_vertices_neighbor(self, vertices_tri, idx_tri, idx_edge):
        """Computes the vertex coordinates of the neighboring triangle across IDX_EDGE given the coordinates of IDX_TRI."""
        vertex_start = vertices_tri[idx_edge]
        vertex_end = vertices_tri[(idx_edge + 1) % 3]

        idx_tri_opp, idx_edge_opp = self.gluings[(idx_tri, idx_edge)]
        tri_opp = self.triangles[idx_tri_opp]
        vertices_tri_opp = {idx_edge_opp: vertex_end,
                            (idx_edge_opp + 1) % 3: vertex_start,
                            (idx_edge_opp + 2) % 3: vertex_start + tri_opp[(idx_edge_opp + 1) % 3]}
        return [vertices_tri_opp[k] for k in range(3)]

    def plot(self):
        # Keys are triangle IDs, and the corresponding values are the vertices of that triangle.
        tris_seen = {0: self.triangles[0].vertices()}
        # Start from the zeroeth triangle in your triangulation.
        tris_to_visit = deque([0])

        plots_tris = self.triangles[0].plot()

        while tris_to_visit:
            idx_tri_curr = tris_to_visit.pop()
            for idx_edge in range(3):
                idx_tri_nbr, _ = self.gluings[(idx_tri_curr, idx_edge)]
                if idx_tri_nbr not in tris_seen:
                    vertices_curr = tris_seen[idx_tri_curr]
                    tris_seen[idx_tri_nbr] = self.get_vertices_neighbor(vertices_curr, idx_tri_curr, idx_edge)
                    plots_tris += self.triangles[idx_tri_nbr].plot(tris_seen[idx_tri_nbr][0])
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

        return plots_tris + plots_labels_tri + plots_labels_edge

    @staticmethod
    def union(tn1, tn2):
        ts = tn1.triangles + tn2.triangles
        tn2_gluings_shifted = {(i1 + len(tn1.triangles), j1): (i2 + len(tn1.triangles), j2)
                               for (i1, j1), (i2, j2) in tn2.gluings.items()}
        gluings = {**tn1.gluings, **tn2_gluings_shifted}
        return Triangulation(ts, gluings)

    @staticmethod
    def merge(tn1, tn2):
        tn = Triangulation.union(tn1, tn2)

        edges_right = {}
        for idx_p, t in enumerate(tn.triangles[len(tn1.triangles):]):
            for idx_e, x in enumerate(t):
                x.set_immutable()
                if x not in tn.gluings:
                    edges_right[x] = (idx_p + len(tn1.triangles), idx_e)

        gluings_new = {}
        for idx_p, t in enumerate(tn.triangles[:len(tn1.triangles)]):
            for idx_e, edge in enumerate(t):
                if (idx_p, idx_e) not in tn.gluings:
                    edge_opp = -edge
                    edge_opp.set_immutable()
                    gluings_new[(idx_p, idx_e)] = edges_right[edge_opp]
                    gluings_new[edges_right[edge_opp]] = (idx_p, idx_e)

        return Triangulation(tn.triangles, {**tn.gluings, **gluings_new})

    @classmethod
    def convex_polygon(cls, edges):
        if len(edges) == 3:
            return Triangulation([Triangle(edges[0], edges[1], edges[2])], gluings={})

        tris = [Triangle(edges[0], edges[1], -(edges[0] + edges[1]))]

        for k in range(2, len(edges) - 2):
            tris.append(Triangle(-tris[-1][2], edges[k], -edges[k] + tris[-1][2]))

        tris.append(Triangle(-(edges[-2] + edges[-1]), edges[-2], edges[-1]))

        gluings = {(k, 2): (k + 1, 0) for k in range(len(edges) - 3)}
        gluings.update({v: k for k, v in gluings.items()})

        return Triangulation(tris, gluings)

    @property
    def area(self):
        return sum(t.area for t in self.triangles)

    @property
    def points_marked(self):
        return {idx: self.triangles[idx].points_marked
                for idx in range(len(self.triangles))}

    @property
    def lines_marked(self):
        return {idx: self.triangles[idx].lines_marked
                for idx in range(len(self.triangles))}


if __name__ == "__main__":
    X = Triangulation.regular_octagon()
    r = X.idr
