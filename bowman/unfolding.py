from sage.all import *

from bowman.triangle import Triangle
from bowman.triangulation import Triangulation


def get_triangle_from_angles(alphas):
    tans = [tan(alphas[0]), tan(alphas[1])]
    _, (a0, a1), phi = number_field_elements_from_algebraics(tans, minimal=True)
    u, v = phi(a0), phi(a1)

    x = v / (u + v)
    y = u * v / (u + v)

    e0 = sage.all.vector([1, 0])
    e2 = -sage.all.vector([x, y])
    return Triangle(e0, -e0 - e2, e2)


def get_genus(a, b, c):
    d = a + b + c
    alpha0, alpha1, alpha2 = QQ(a / d), QQ(b / d), QQ(c / d)
    return 1 + d / 2 * (1 - sum(1 / q.denominator() for q in [alpha0, alpha1, alpha2]))


def enumerate_triangles(height, num_fake=None):
    for (a, b, c) in [(a, b, c)
                      for a in range(1, height)
                      for b in range(a, height)
                      for c in range(b, height)
                      if a <= b <= c
                      if a + b + c <= height
                      if gcd(a, gcd(b, c)) == 1
                      if get_genus(a, b, c) >= 2]:
        d = a + b + c
        q0, q1, q2 = QQ(a / d), QQ(b / d), QQ(c / d)
        data_cone_pts = [f"{d // q.denominator()} points of multiplicity {q.numerator()}" for q in [q0, q1, q2]]
        data_total = f"(q0, q1, q2)={(q0, q1, q2)}, {data_cone_pts}"
        num_numerators_bad = sum(q.numerator() == 1 for q in [q0, q1, q2])

        if num_fake is None:
            yield (a, b, c), data_total
        elif num_fake == num_numerators_bad:
            yield (a, b, c), data_total


def rotation_matrix(theta):
    c = QQbar(cos(theta))
    s = QQbar(sin(theta))
    return sage.all.matrix([[c, -s], [s, c]])


def reflection_matrix(w):
    """Get reflection matrix across Span(w)"""
    w_perp = sage.all.matrix([[0, -1], [1, 0]]) * w
    change_basis = sage.all.matrix([w, w_perp]).transpose()  # COB: (e0, e1) --> (v, v_perp)
    return change_basis * sage.all.matrix([[1, 0], [0, -1]]) * change_basis.inverse()


def _develop_cone_point(triangle, alpha):
    """Rotate triangle about the given vertex and return orbit"""
    p, q = QQ(alpha / pi).numerator(), QQ(alpha / pi).denominator()
    m = rotation_matrix(alpha)
    if (2 * q) % p == 0:
        edge_succ = triangle[1]
        num_sides = QQ(p / (2 * q)).denominator()
        edges = []
        m1 = sage.all.identity_matrix(2)
        for _ in range(num_sides):
            edges.append(m1 * edge_succ)
            m1 *= m
        return Triangulation.convex_polygon(edges)
    else:
        triangles = []
        num_triangles = 2 * q if p % 2 == 1 else q
        m1 = sage.all.identity_matrix(2)
        for _ in range(num_triangles):
            triangles.append(Triangle(*[m1 * x for x in triangle]))
            m1 *= m
        gluings = {(k, 2): ((k + 1) % num_triangles, 0) for k in range(num_triangles)}
        gluings.update({v: k for k, v in gluings.items()})
        return Triangulation(triangles, gluings)


def _develop_cone_points(tri, alpha, num_copies, rot_mat):
    tns = []
    m = sage.all.identity_matrix(2)
    for k in range(num_copies):
        tri_curr = Triangle(m * tri[0], m * tri[1], m * tri[2])
        tns.append(_develop_cone_point(tri_curr, alpha))
        m *= rot_mat

    tn = Triangulation()
    for tn1 in tns:
        tn = Triangulation.union(tn, tn1)

    return tn


def _triangulate_unfolding(a, b, c):
    """Decompose unfolding of triangle into polygons+gluings with real vertices"""
    d = a + b + c
    alphas = [QQ(a / d) * pi, QQ(b / d) * pi, QQ(c / d) * pi]
    tri = get_triangle_from_angles(alphas)

    refl = sage.all.matrix([[1, 0], [0, -1]])
    e0, e1, e2 = tri
    tri_left = Triangle(-(refl * e2), (refl * e2) - e2, e2)
    tri_right = Triangle(e1, -e1 + (refl * e1), -(refl * e1))

    (p0, q0), (p1, q1), (p2, q2) = [(QQ(x / pi).numerator(), QQ(x / pi).denominator()) for x in alphas]
    num_quads = sage.all.lcm(q0, sage.all.lcm(q1, q2))

    rot_l = rotation_matrix(2 * alphas[0])
    rot_r = rotation_matrix(2 * alphas[1])

    tn_left = _develop_cone_points(tri_left, 2 * alphas[0], num_quads // q0, rot_r)
    tn_right = _develop_cone_points(tri_right, 2 * alphas[1], num_quads // q1, rot_l)

    return Triangulation.merge(tn_left, tn_right).make_nontrivial()


def cos_angle(v, w):
    return v.dot_product(w) / (v.norm() * w.norm())


def _build_field_gothic(c):
    d = 4 * (c ** 2) + 1
    x = PolynomialRing(QQ, 'x').gen()
    poly = (x ** 4 - 2 * (d + 3) * (x ** 2) + (d - 3) ** 2)
    K = NumberField(poly, 'u', embedding=QQ(1))
    u = K.gen()  # u = rtd - rt3
    rt3d = (u ** 2 - 3 - d) / (-2)
    rtd = ((rt3d + d) * u) / (d - 3)
    rt3 = -u + rtd

    return K, [rt3, rtd]


def triangulate_gothic1128(c):
    K, [rt3, rtd] = _build_field_gothic(c)
    y = K(1 / 2) * (rtd - K(2 * c + 1))
    assert y > 0 and not y.is_rational() and y**2 + (2*c + 1)*y + c == 0

    vertices = [sage.all.vector([K(0), K(0)]),
                sage.all.vector([y / 2, -rt3 / 2 * y]),
                sage.all.vector([QQ(3 / 2) * y + QQ(1 / 2), rt3 / 2 * y + rt3 / 2]),
                sage.all.vector([K(-1), K(0)])]

    [e0, e1, e2, e3] = [vertices[1] - vertices[0], vertices[2] - vertices[1],
                        vertices[3] - vertices[2], vertices[0] - vertices[3]]

    e0_mirror, e3_mirror = reflection_matrix(e1) * -e0, reflection_matrix(-e2) * e3

    t0, t1 = Triangle(e3, -e3 + e3_mirror, -e3_mirror), Triangle(e0_mirror, -e0_mirror - e0, e0)

    rot = sage.all.matrix([[K(1/2), -rt3/2], [rt3/2, K(1/2)]])  # rotation by pi/3
    tn0 = Triangulation.convex_polygon([(rot ** k) * t0[1] for k in range(6)])
    tn1 = Triangulation.union(Triangulation.convex_polygon([rot ** (2 * k) * t1[1] for k in range(3)]),
                              Triangulation.convex_polygon([rot ** (2 * k) * (rot * t1[1]) for k in range(3)]))
    tn2 = Triangulation.convex_polygon([x for k in range(6)
                                        for x in [(rot ** k) * -t0[1], (rot ** k) * -t1[1]]])

    tn = Triangulation.merge(tn1, tn2)
    tn = Triangulation.merge(tn0, tn)

    return tn.make_nontrivial()


if __name__ == "__main__":
    gothic1228 = triangulate_gothic1128(QQ(-1))
    fund_dom = gothic1228.generators_veech
