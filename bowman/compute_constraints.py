import itertools

from sage.all import *


def compute_constraints_a_plus(d):
    assert d >= 17 and not ZZ(d).is_square(), "disc must be >= 17 and not a square"
    k = QuadraticField(d)
    rt_d = k.gen()
    if d % 8 == 0:
        w, lmbd = k(d) / 8, rt_d / 2
    elif d % 8 == 1:
        w, lmbd = k(d - 1) / 8, (-1 + rt_d) / 2
    elif d % 8 == 4:
        w, lmbd = k(d - 4) / 8, (-2 + rt_d) / 2
    else:
        assert False, "disc must be 0, 1, or 4 mod 8"

    for ell, height_v_new in itertools.product([k(0), k(1)], [w - lmbd, lmbd]):
        constraints_h0_v0_horizontal = sage.all.matrix(k, [[0, 1, 0],
                                                           [1 / (w - lmbd), 0, 0],
                                                           [1 / (w - lmbd), w / (w - lmbd), (-ell * w) / height_v_new]])
        try:
            [[_, a01, _], [a10, _, _], [a20, a21, a22]] = (constraints_h0_v0_horizontal ** -1)
        except ZeroDivisionError:
            print("broke on", ell, height_v_new)
        else:
            print(a22[1])


if __name__ == "__main__":
    compute_constraints_a_plus(24)
