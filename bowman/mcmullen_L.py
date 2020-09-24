from sage.all import *

from bowman.triangulation import Triangle, Triangulation


def gen_mcmullen_l(a, b):
    def triangulate_rectangle(base, height):
        triangle_lower = Triangle(sage.all.vector([0, height]),
                                  sage.all.vector([-base, -height]),
                                  sage.all.vector([base, 0]))
        triangle_upper = Triangle(sage.all.vector([0, -height]),
                                  sage.all.vector([base, height]),
                                  sage.all.vector([-base, 0]))
        return [triangle_lower, triangle_upper]

    if not (a > 1 and b > 1):
        raise ValueError("Need to have a, b > 1")
    elif a.parent() != b.parent():
        raise ValueError("a, b need to come from same field")
    triangles = [*triangulate_rectangle(1, 1),
                 *triangulate_rectangle(b - 1, 1),
                 *triangulate_rectangle(1, a - 1)]

    gluings = {(0, 0): (3, 0), (0, 1): (1, 1), (0, 2): (5, 2),
               (1, 0): (2, 0), (1, 2): (4, 2), (2, 1): (3, 1),
               (2, 2): (3, 2), (4, 0): (5, 0), (4, 1): (5, 1)}
    gluings.update({v: k for k, v in gluings.items()})

    return Triangulation(triangles, gluings, a.parent())


if __name__ == "__main__":
    X = gen_mcmullen_l(QQ(5 / 4), QQ(17 / 10))
    print(Triangulation.mcmullen_l(QQ(5 / 4), QQ(17 / 10)))
