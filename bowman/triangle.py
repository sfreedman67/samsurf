from collections import namedtuple
import functools

from sage.all import *


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
        return sage.all.polygon2d(self.vertices(basepoint), fill=False).plot()

    def vertices(self, basepoint=sage.all.zero_vector(2)):
        return [basepoint, basepoint + self.v0, basepoint - self.v2]

    def apply_matrix(self, m):
        w0 = m * self.v0
        w1 = m * self.v1
        w2 = -(w0 + w1)
        return Triangle(w0, w1, w2)

    def __hash__(self):
        return hash(tuple(coord for vertex in self.vertices() for coord in vertex))

    @property
    def area(self):
        return QQ(1/2) * abs(sage.all.matrix([self.v0, -self.v2]).determinant())
