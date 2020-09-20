from sage.all import *

from collections import namedtuple

from bowman.triangulation import Triangulation

PointR2 = namedtuple("PointR2", ["x", "y"])


def translate(p, v):
    return PointR2(p.x + v[0], p.y + v[1])


TriangleRendered = namedtuple("TriangleRendered", ["p0", "p1", "p2"])


def plot_tris_rendered(ts_r):
    return sum(sage.plot.polygon.polygon(t_r, fill=False) for t_r in ts_r)


def render_tri_on_edge(t, initial, final):
    for idx, v_nbr in enumerate(t):
        if v_nbr == -sage.all.vector([final.x - initial.x, final.y - initial.y]):
            vert_new = translate(initial, t[(idx + 1) % 3])
            verts_nbr = {idx: final,
                         (idx + 1) % 3: initial,
                         (idx + 2) % 3: vert_new}
            return TriangleRendered(verts_nbr[0], verts_nbr[1], verts_nbr[2])


def render_tri_at_origin(t):
    origin = PointR2(0, 0)
    p1 = translate(origin, t[0])
    p2 = translate(p1, t[1])
    return TriangleRendered(origin, p1, p2)


def render_across_edge(trin, edge_comb, t_r):
    _, idx_edge = edge_comb
    idx_tri_opp, _ = trin.gluings[edge_comb]
    return render_tri_on_edge(X.triangles[idx_tri_opp],
                              t_r[idx_edge],
                              t_r[(idx_edge + 1) % 3])


if __name__ == "__main__":
    # Let's produce a rendering for the regular octagon
    X = Triangulation.regular_octagon()
    # First, render t0
    t0_r = render_tri_at_origin(X.triangles[0])
    # Next, render its neighbors
    t5_r, t2_r, t3_r = (render_across_edge(X, (0, k), t0_r) for k in range(3))
    # Go to t5 and render its neighbors
    _, t4_r, t1_r = (render_across_edge(X, (5, k), t5_r) for k in range(3))
    # And that's it! We can automate with a tree traversal
    # But what about
    plot_tris_rendered([t0_r, t5_r, t2_r, t3_r, t4_r, t1_r]).show()
