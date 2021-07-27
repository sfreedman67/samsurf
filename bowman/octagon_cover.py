import itertools
from collections import namedtuple

from bowman.triangle import Triangle
from bowman.triangulation import Triangulation

from sage.all import *

class OctagonCover(Triangulation):

    def __make_octagon_cover(self, triangulation):
        #self = OctagonCover(degree, (0, 0, 0, 0, 0, 0, 0, 0, 0))
        self.triangles = triangulation.triangles
        self.gluings = triangulation.gluings
        return self

    def __init__(self, degree, generator_images=None):

        def holonomy_representation(generator_images):
            """Constructs the tuple to be fed into the OctagonCover constructor.

            generator_images := a tuple of nine nonnegative integers. Of the first
                                eight, the evens correspond to the four standard
                                generators of the fundamental group of the unpunctured
                                surface, and the odds correspond to those additional
                                generators adjoined as a result of the four punctures
                                on the edges. Finally, the last component indicates the
                                holonomy of the final puncture in the center of the
                                octagon."""
            holonomy_rep = tuple(generator_images[:-1])
            for i in range(len(holonomy_rep) // 2):
                holonomy_rep = holonomy_rep + (-holonomy_rep[2 * i + 1], -holonomy_rep[2 * i])
            holonomy_rep = holonomy_rep + (generator_images[-1], -generator_images[-1])
            return holonomy_rep

        self.degree = degree
        holonomy_rep = holonomy_representation(generator_images)

        # Step 1: Construct the triangles.
        k = sage.all.QuadraticField(2)
        sqrt2 = k.gen()

        h_long = sage.all.vector([Rational('1/2'), 0])
        h_short = sage.all.vector([Rational('1/4') * sqrt2, 0])
        v_long = sage.all.vector([0, Rational('1/2')])
        v_short = sage.all.vector([0, Rational('1/4') * sqrt2])
        ne_diagonal = h_long + h_short + v_long + v_short
        nw_diagonal = -(h_long + h_short) + v_long + v_short

        triangles_one_octagon = (Triangle(h_long + 2 * h_short + v_long, -h_short + v_short, -ne_diagonal),
                                 Triangle(ne_diagonal, -h_short + v_short, -(h_long + v_long + 2 * v_short)),
                                 Triangle(h_long + v_long + 2 * v_short, -h_long, -(v_long + 2 * v_short)),
                                 Triangle(v_long + 2 * v_short, -h_long, h_long - (v_long + 2 * v_short)),
                                 Triangle(-h_long + v_long + 2 * v_short, -h_short - v_short, -nw_diagonal),
                                 Triangle(nw_diagonal, -h_short - v_short, h_long + 2 * h_short - v_long),
                                 Triangle(-(h_long + 2 * h_short) + v_long, -v_long, h_long + 2 * h_short),
                                 Triangle(-(h_long + 2 * h_short), -v_long, h_long + 2 * h_short + v_long),
                                 Triangle(-(h_long + 2 * h_short) - v_long, h_short - v_short, ne_diagonal),
                                 Triangle(-ne_diagonal, h_short - v_short, h_long + v_long + 2 * v_short),
                                 Triangle(-h_long - (v_long + 2 * v_short), h_long, v_long + 2 * v_short),
                                 Triangle(-(v_long + 2 * v_short), h_long, -h_long + v_long + 2 * v_short),
                                 Triangle(h_long - (v_long + 2 * v_short), h_short + v_short, nw_diagonal),
                                 Triangle(-nw_diagonal, h_short + v_short, -(h_long + 2 * h_short) + v_long),
                                 Triangle(h_long + 2 * h_short - v_long, v_long, -(h_long + 2 * h_short)),
                                 Triangle(h_long + 2 * h_short, v_long, -(h_long + 2 * h_short) - v_long))

        triangles = tuple()
        for i in range(degree):
            triangles = triangles + triangles_one_octagon

        # Step 2: Define edge pairings on the base.
        gluings_outer = dict()

        for i in range(4):
            gluings_outer[2 * i] = 2 * (i + 4) + 1
            gluings_outer[2 * i + 1] = 2 * (i + 4)

        gluings_outer.update({v: k for k, v in gluings_outer.items()})

        # Step 3: Figure out the gluings in the covering space.
        gluings_cover = dict()

        for i in range(degree):
            # a. Piece together the interior edges.
            for j in range(16 - 1):
                gluings_cover[(16 * i + j, 2)] = (16 * i + j + 1, 0)
            # b. Take care of the interior slit.
            gluings_cover[(16 * i + 15, 2)] = (16 * ((i + holonomy_rep[16]) % degree), 0)
            # c. Piece together the outer edges.
            for j in range(16):
                tri_in = 16 * ((i + holonomy_rep[j]) % degree) + gluings_outer[j]
                gluings_cover[(16 * i + j, 1)] = (tri_in, 1)

        gluings_cover.update({v: k for k, v in gluings_cover.items()})

        super().__init__(triangles, gluings_cover)

    def flow_permutation(self, base_tri_id, base_coords, direction):
        """Given a point in the base and a cylinder direction, computes the
        action of the flow on the fiber above the point."""
        perm = dict()
        for i in range(self.degree):
            is_found_preimage = False
            tri_id, coords, _ = self.step_flow(16 * i + base_tri_id, base_coords, direction)
            while not is_found_preimage:
                # 1. Flow forward until we arrive at a congruent triangle.
                while tri_id % 16 != base_tri_id:
                    tri_id, coords, _ = self.step_flow(tri_id, coords, direction)

                # 2. Determine whether we run into a point in the fiber.
                tri = self.triangles[tri_id]
                start_pos = coords[1] * tri[0] - coords[2] * tri[2]
                end_pos = base_coords[1] * tri[0] - base_coords[2] * tri[2]
                displacement = end_pos - start_pos
                nonzero_component = 0
                if direction[nonzero_component].is_zero():
                    nonzero_component = 1
                ratio = displacement[nonzero_component] / direction[nonzero_component]
                if displacement == ratio * direction:
                    is_found_preimage = True
            perm[i] = tri_id // 16
        return perm

    def mark_line(self, triangle_id, start_coords, end_coords, rgbcolor):
        return self.__make_octagon_cover(super().mark_line(triangle_id, start_coords, end_coords, rgbcolor))

    def mark_point(self, triangle_id, coords, rgbcolor):
        return self.__make_octagon_cover(super().mark_point(triangle_id, coords, rgbcolor))

    def mark_flow(self, start_tri_id, start_coords, direction, rgbcolor):
        return self.__make_octagon_cover(super().mark_flow(start_tri_id, start_coords, direction, rgbcolor))

    def plot(self):
        """Plot as individual octagons"""

        tris_seen = {0: self.triangles[0].vertices()}
        plots_tris_temp = []

        #for each octagon
        for i in range(self.degree):
            #produce an empty graphic to hold the octagon
            octagon_curr = Graphics()
            for j in range(16):
                index = i * 16 + j
                #add triangle at each index
                tris_seen[index] = self.triangles[index].vertices()
                octagon_curr += self.triangles[index].plot()
            #add octagon to plot
            plots_tris_temp.append(octagon_curr)

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

        for i in range(self.degree):
            for j in range(16):
                idx = i*16 + j
                (a, b, c) = tris_seen[idx]
                plots_tris_temp[i] += sage.all.text(str(idx), center(tris_seen[idx]), fontsize=14, color='orange').plot()
                plots_tris_temp[i] += sum(sage.all.text(str(idx), midpoint(v1, v2) + displacement(v1, v2)).plot()
                                         for idx, (v1, v2) in [(0, (a, b)), (1, (b, c)), (2, (c, a))])

        return graphics_array(plots_tris_temp, 1, self.degree)