from dataclasses import dataclass, field

from sage.plot.polygon import polygon2d
from sage.plot.line import line2d
from sage.plot.point import point2d
from sage.rings.number_field.number_field import QuadraticField
from sage.modules.free_module_element import vector


@dataclass(frozen=True)
class Rectangle:
    field: ...
    base: ...
    height: ...

    @property
    def coords(self):
        return [vector(self.field, [0, 0]),
                vector(self.field, [self.base, 0]),
                vector(self.field, [self.base, self.height]),
                vector(self.field, [0, self.height])]


@dataclass(frozen=True)
class TV():
    field: ...
    rectangles: tuple[Rectangle] = field(default_factory=tuple)
    right: tuple[int] = field(default_factory=tuple)
    top: tuple[int] = field(default_factory=tuple)
    instrs: ... = None

    @classmethod
    def prym2(cls, w, e):
        d = e**2 + 4 * w * 1
        assert d >= 5

        k = QuadraticField(d)
        rtd = k.gen()
        lmbd = (e + rtd) / 2

        sq = Rectangle(k, lmbd, lmbd)
        r_short = Rectangle(k, lmbd, k(1))
        r_long = Rectangle(k, k(w) - lmbd, k(1))

        right = (1, 0, 2)
        top = (2, 1, 0)

        instrs = [(0, 'r', 1), (0, 'u', 2)]

        return TV(k, [r_short, r_long, sq], right, top, instrs)

    @classmethod
    def prym4(cls, w, e):
        d = e**2 + 4 * w
        assert (e + 4)**2 < d
        k = QuadraticField(d)
        rtd = k.gen()
        lmbd = (e + rtd)/2

        sq = Rectangle(k, lmbd/2, lmbd/2)
        r_mid = Rectangle(k, lmbd/2, k(1/2))
        r_thin = Rectangle(k, k(w/2) - lmbd, k(1/2))

        right = (0, 2, 3, 1, 5, 6, 4, 7)
        top = (3, 1, 4, 5, 7, 0, 6, 2)

        instrs = [(0, 'u', 3), (3, 'l', 2), (2, 'l', 1), (3, 'u', 5),
                  (5, 'r', 6), (5, 'l', 4), (4, 'u', 7)]

        return TV(k, [sq, r_thin, r_mid, r_mid, r_mid, r_mid, r_thin, sq],
                  right, top, instrs)

    def _to_samsurf(self):
        from context import samsurf
        from triangle import Triangle
        from triangulation import Triangulation
        tris = []
        e1 = vector(self.field, [1, 0])
        e2 = vector(self.field, [0, 1]) 
        for r in self.rectangles:
            t1 = Triangle(r.base * e1,
                          r.height * e2,
                          -r.base * e1 + -r.height * e2)
            t2 = Triangle(-r.base * e1,
                          -r.height * e2,
                          r.base * e1 + r.height * e2)
            tris.extend([t1, t2])

        glues = {}
        for idx_rect, r in enumerate(self.rectangles):
            glues[(2 * idx_rect, 0)] = (2 * self.bottom[idx_rect] + 1, 0)
            glues[(2 * idx_rect, 1)] = (2 * self.right[idx_rect] + 1, 1)
            glues[(2 * idx_rect, 2)] = (2 * idx_rect + 1, 2)

        glues.update({v: k for k, v in glues.items()})
        return Triangulation(tris, glues)


    @property
    def _corners(self):
        if self.instrs is None:
            raise ValueError

        corners = [vector(self.field, [0, 0])
                   for _ in range(len(self.rectangles))]
        for idx0, dirn, idx1 in self.instrs:
            r0, r1 = self.rectangles[idx0], self.rectangles[idx1]
            translation = {'r': vector(self.field, [r0.base, 0]),
                           'l': vector(self.field, [-r1.base, 0]),
                           'u': vector(self.field, [0, r0.height]),
                           'd': vector(self.field, [0, -r1.height])}
            corners[idx1] = corners[idx0] + translation[dirn]
        return corners

    def plot(self, points=None, segments=None):
        if self.instrs is None:
            raise ValueError

        rects = []
        for idx, corner in enumerate(self._corners):
            coords1 = [corner + coord
                       for coord in self.rectangles[idx].coords]
            rects.append(polygon2d(coords1, fill=False, color='black'))
        plt = sum(rects)

        if points is not None:
            for idx_rect, point, c in points:
                corner = self._corners[idx_rect]
                plt += point2d([corner + point], color=c, size=18)

        if segments is not None:
            for idx_rect, seg, c in segments:
                corner = self._corners[idx_rect]
                plt += line2d([corner + seg[0], corner + seg[1]],
                              color=c)
        return plt

    def decompose(self, m, color="black"):
        return [x for idx in range(len(self.rectangles))
                for x in self.flow(idx, vector(self.field, [0, 0]), m, color)]

    def flow(self, idx_init, pt_init, m, color="black"):
        data = []
        idx, pt0 = idx_init, pt_init
        while True:
            rect = self.rectangles[idx]
            pt1 = TV.step_flow(rect, pt0, m)
            x1, y1 = pt1
            if idx == idx_init and TV._is_between_strict(pt0, pt_init, pt1):
                data.append((idx, (pt0, pt_init), color))
                return data

            data.append((idx, (pt0, pt1), color))

            if (x1, y1) == (rect.base, rect.height):
                if self.top[self.right[idx]] != self.right[self.top[idx]]:
                    return data
                idx = self.top[self.right[idx]]
                pt0 = vector(self.field, [0, 0])
            elif x1 == rect.base:
                idx = self.right[idx]
                pt0 = vector(self.field, [0, y1])
            elif y1 == rect.height:
                idx = self.top[idx]
                pt0 = vector(self.field, [x1, 0])
            else:
                assert False

            if idx == idx_init and pt0 == pt_init:
                return data

    @staticmethod
    def step_flow(rect, pt, m):
        assert m >= 0
        x, y = pt
        m_corner = (rect.height - y)/(rect.base - x)
        if m < m_corner:
            # intersects right side
            return vector(rect.field, [rect.base, y + m * (rect.base - x)])
        elif m == m_corner:
            # intersects corner
            return vector(rect.field, [rect.base, rect.height])
        else:
            # intersects top
            return vector(rect.field, [x + (rect.height - y)/m, rect.height])

    @staticmethod
    def _is_between_strict(p0, p1, p2):
        d01 = (p1 - p0).norm()
        d12 = (p2 - p1).norm()
        d02 = (p2 - p0).norm()

        return all([d01 != 0, d12 != 0, d01 + d12 == d02])
