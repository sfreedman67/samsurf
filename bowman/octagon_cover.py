import itertools
from collections import namedtuple

from bowman.triangle import Triangle
from bowman.triangulation import Triangulation

from sage.all import *

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

class OctagonCover(Triangulation):

    def __init__(self, degree, holonomy_rep):
        self.degree = degree

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

    def plot(self):
        return super().octagon_plot(self.degree)