import itertools

from sage.all import *

from bowman.triangle import Triangle
from bowman.triangulation import Triangulation


def get_triangle_from_angles(alphas):
    tans = [tan(alphas[0]), tan(alphas[1])]
    K, (tan0, tan1), _ = number_field_elements_from_algebraics(tans, minimal=True, embedded=True)

    e0 = sage.all.vector([tan0 + tan1, 0])
    e2 = -sage.all.vector([tan1, tan0 * tan1])
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
    trig = [cos(theta), sin(theta)]
    _, (c, s), _ = number_field_elements_from_algebraics(trig, minimal=True, embedded=True)
    return sage.all.matrix([[c, -s], [s, c]])


def _develop_cone_point(triangle, alpha, m):
    """Rotate triangle about the given vertex and return orbit"""
    p, q = QQ(alpha / pi).numerator(), QQ(alpha / pi).denominator()
    # m = rotation_matrix(alpha)
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


def _develop_cone_points(tri, alpha, rot_vertex, num_copies, rot_mat):
    tns = []
    m = sage.all.identity_matrix(2)
    for k in range(num_copies):
        tri_curr = Triangle(m * tri[0], m * tri[1], m * tri[2])
        tns.append(_develop_cone_point(tri_curr, alpha, rot_vertex))
        m *= rot_mat

    tn = Triangulation()
    for tn1 in tns:
        tn = Triangulation.union(tn, tn1)

    return tn


def _triangulate_unfolding(a, b, c):
    """Decompose unfolding of triangle into polygons+gluings with real vertices"""
    d = a + b + c
    alphas = [QQ(a / d) * pi, QQ(b / d) * pi, QQ(c / d) * pi]
    trig = [cos(alphas[0]), sin(alphas[0]), cos(alphas[1]), sin(alphas[1])]
    _, (c0, s0, c1, s1), _ = number_field_elements_from_algebraics(trig, minimal=True, embedded=True)
    tan0, tan1 = s0 / c0, s1 / c1

    e0 = sage.all.vector([tan0 + tan1, 0])
    e2 = -sage.all.vector([tan1, tan0 * tan1])
    tri = Triangle(e0, -e0 - e2, e2)

    refl = sage.all.matrix([[QQ(1), QQ(0)], [QQ(0), QQ(-1)]])
    e0, e1, e2 = tri
    tri_left = Triangle(-(refl * e2), (refl * e2) - e2, e2)
    tri_right = Triangle(e1, -e1 + (refl * e1), -(refl * e1))

    (p0, q0), (p1, q1), (p2, q2) = [(QQ(x / pi).numerator(), QQ(x / pi).denominator()) for x in alphas]
    num_quads = sage.all.lcm(q0, sage.all.lcm(q1, q2))

    two_s0, two_s1 = 2 * s0 * c0, 2 * s1 * c1
    two_c0, two_c1 = c0 ** 2 - s0 ** 2, c1 ** 2 - s1 ** 2
    rot_l = sage.all.matrix([[two_c0, -two_s0], [two_s0, two_c0]])
    rot_r = sage.all.matrix([[two_c1, -two_s1], [two_s1, two_c1]])

    tn_left = _develop_cone_points(tri_left, QQ(2) * alphas[0], rot_l, num_quads // q0, rot_r)
    tn_right = _develop_cone_points(tri_right, QQ(2) * alphas[1], rot_r, num_quads // q1, rot_l)

    return Triangulation.merge(tn_left, tn_right).make_nontrivial()


def _build_field_gothic1128(c):
    d = 4 * (c ** 2) + 1
    x = PolynomialRing(QQ, 'x').gen()
    poly = (x ** 4 - 2 * (d + 3) * (x ** 2) + (d - 3) ** 2)
    K = NumberField(poly, 'u', embedding=QQ(1))
    u = K.gen()  # u = rtd - rt3
    rt3d = (u ** 2 - 3 - d) / (-2)
    rtd = ((rt3d + d) * u) / (d - 3)
    rt3 = -u + rtd

    return K, [rt3 if rt3 > 0 else -rt3, rtd if rtd > 0 else -rtd]


def _build_field_gothic(c, angles):
    if angles == "1119":
        d = 9 * c ** 2 + 2 * c + 1
    elif angles == "1128":
        d = 4 * (c ** 2) + 1
    else:
        raise ValueError("Mode is either \"1128\" or \"1119\"")
    x = PolynomialRing(QQ, 'x').gen()
    poly = (x ** 4 - 2 * (d + 3) * (x ** 2) + (d - 3) ** 2)
    K = NumberField(poly, 'u', embedding=QQ(1))
    u = K.gen()  # u = rtd - rt3
    rt3d = (u ** 2 - 3 - d) / (-2)
    rtd = ((rt3d + d) * u) / (d - 3)
    rt3 = -u + rtd

    return K, [rt3 if rt3 > 0 else -rt3, rtd if rtd > 0 else -rtd]


def triangulate_gothic1128(c):
    if c >= 0:
        raise ValueError("c must be negative")
    elif (4 * c ** 2 + 1).is_square():
        raise ValueError("y must be irrational")

    K, [rt3, rtd] = _build_field_gothic(c, angles="1128")
    y = - QQ(c) - QQ(1 / 2) + rtd / 2

    vertices = [sage.all.vector([QQ(0), QQ(0)]),
                sage.all.vector([y / 2, -rt3 / 2 * y]),
                sage.all.vector([QQ(3 / 2) * y + QQ(1 / 2), rt3 / 2 * y + rt3 / 2]),
                sage.all.vector([QQ(-1), QQ(0)])]

    e0, e3 = vertices[1] - vertices[0], vertices[0] - vertices[3]

    rot = sage.all.matrix([[QQ(1 / 2), -rt3 / 2], [rt3 / 2, QQ(1 / 2)]])  # rotation by pi/3
    e0_mirror, e3_mirror = rot ** (-2) * -e0, rot * e3
    t0, t1 = Triangle(e3, -e3 + e3_mirror, -e3_mirror), Triangle(e0_mirror, -e0_mirror - e0, e0)

    tn0 = Triangulation.convex_polygon([(rot ** k) * t0[1] for k in range(6)])
    tn1 = Triangulation.union(Triangulation.convex_polygon([rot ** (2 * k) * t1[1] for k in range(3)]),
                              Triangulation.convex_polygon([rot ** (2 * k) * (rot * t1[1]) for k in range(3)]))
    tn2 = Triangulation.convex_polygon([x for k in range(6)
                                        for x in [(rot ** k) * -t0[1], (rot ** k) * -t1[1]]])
    tn = Triangulation.merge(tn1, tn2)
    tn = Triangulation.merge(tn0, tn)
    return tn.make_nontrivial()


def triangulate_gothic1119(c):
    if c >= 0:
        raise ValueError("c must be negative")
    elif (9 * c ** 2 + 2 * c + 1).is_square():
        raise ValueError("y must be irrational")

    K, [rt3, rtd] = _build_field_gothic(c, angles="1119")
    y = (-QQ(3) * c - QQ(1) + rtd) / 2

    vertices = [sage.all.vector([QQ(0), QQ(0)]),
                sage.all.vector([QQ(0), -rt3 * y]),
                sage.all.vector([QQ(3 / 2) * y + QQ(1 / 2), (QQ(1 / 2) * y + QQ(1 / 2)) * rt3]),
                sage.all.vector([QQ(-1), QQ(0)])]

    e0, e3 = vertices[1] - vertices[0], vertices[0] - vertices[3]

    rot = sage.all.matrix([[QQ(1 / 2), -rt3 / 2], [rt3 / 2, QQ(1 / 2)]])  # rotation by pi/3
    e0_mirror, e3_mirror = rot ** (-1) * -e0, rot * e3
    t0, t1 = Triangle(e3, -e3 + e3_mirror, -e3_mirror), Triangle(e0_mirror, -e0_mirror - e0, e0)

    tn0 = Triangulation.convex_polygon([(rot ** k) * t0[1] for k in range(6)])
    tn1 = Triangulation.convex_polygon([(rot ** k) * t1[1] for k in range(6)])
    tn2 = Triangulation.convex_polygon([x for k in range(6)
                                        for x in [(rot ** k) * -t0[1], (rot ** k) * -t1[1]]])
    tn = Triangulation.merge(tn1, tn2)
    tn = Triangulation.merge(tn0, tn)
    return tn.make_nontrivial()


def get_chis(max_height):
    for height in range(1, max_height + 1):
        for n1, n2 in [(a, height - a) for a in range(1, height) if gcd(a, height - a) == 1]:
            c = -QQ(n1 / n2)
            gothic1228 = triangulate_gothic1128(c)
            fund_dom = gothic1228.generators_veech
            print(f"c={c}, num_cusps={fund_dom.cusps},"
                  f"chi_orb={fund_dom.chi_orb}, orbifold_points={fund_dom.points_orbifold}")


def get_proportions(max_height):
    from bowman.algo import BoundaryIsNotZippedPairwiseError
    for rat in [QQ(-1 / t) for t in range(1, max_height)]:
        print(rat)
        X = triangulate_gothic1128(rat)
        try:
            f = X.generators_veech
        except BoundaryIsNotZippedPairwiseError:
            print(f"rat={rat} is bad")
        else:
            for code, pct in sorted(f.proportions.items(), key=lambda x: -x[1]):
                print(code, "{:.5%}".format(float(pct)))
    print("---")


if __name__ == "__main__":
    rats = set()
    depth = 20
    for p, q in itertools.product(range(1, depth), repeat=2):
        r = QQ(-p/q)
        D = ZZ(4 * p**2 + q**2).squarefree_part()
        if D != 1:
            rats.add((r, D))
    for D in sorted(set(disc for _, disc in rats)):
        r = min((rat for rat, disc in rats if disc == D), key=lambda x: max(abs(x.numerator()), abs(x.denominator())))
        print(r, D)
        X = triangulate_gothic1128(r)
        f = X.generators_veech
