from dataclasses import dataclass, field

import numpy as np

from sage.all import vector
from sage.plot.polygon import polygon2d
from sage.plot.line import line2d
from sage.plot.point import point2d
from sage.rings.number_field.number_field import QuadraticField

@dataclass(frozen=True)
class Rectangle:
    base: ...
    height: ...
        
    @property
    def coords(self):
        return [np.zeros(2),
                np.array([self.base, 0]),
                np.array([self.base, self.height]),
                np.array([0, self.height])]

@dataclass(frozen=True)
class TV():
    field: ...
    rectangles: tuple[Rectangle] = field(default_factory=tuple)
    right: tuple[int] = field(default_factory=tuple)
    top: tuple[int] = field(default_factory=tuple)


    @classmethod
    def prym2(cls, d):
        k = QuadraticField(d)
        rtd = k.gen()
        if d % 4 == 0:
            w, e = d / 4, 0
        elif d % 4 == 1:
            w, e = (d - 1) / 4, -1
        else:
            assert False
        lmbd = (e + rtd) / 2

        sq = Rectangle(lmbd, lmbd)
        r_short = Rectangle(lmbd, k(1))
        r_long = Rectangle(w - lmbd, k(1))

        right = (1, 0, 2)
        top = (2, 1, 0)

        return TV(k, [r_short, r_long, sq], right, top)

    @classmethod
    def prym4(cls, d):
        k = QuadraticField(d)
        rtd = k.gen()
        if d % 4 == 0:
            w, e = d // 4, 0
        elif d % 4 == 1:
            w, e = (d - 1) // 4, -1
        else:
            assert False
        lmbd = (e + rtd)/2

        sq = Rectangle(lmbd/2, lmbd/2)
        r_mid = Rectangle(lmbd/2, k(1/2))
        r_thin = Rectangle(w/2 - lmbd, k(1/2))

        right = (0, 2, 3, 1, 5, 6, 4, 7)
        top = (3, 1, 4, 5, 7, 0, 6, 2)

        return TV(k, [sq, r_thin, r_mid, r_mid, r_mid, r_mid, r_thin, sq],
                  right, top)

    @property
    def left(self):
        return tuple(self.right.index(k) for k in range(len(self.rectangles)))

    @property
    def bottom(self):
        return tuple(self.top.index(k) for k in range(len(self.rectangles)))

    def _to_samsurf(self):
        from context import bowman
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


    def plot(self, instrs, points=None, segs=None):
        corners = self._get_corners(instrs)

        rects = []
        for idx, corner in enumerate(corners):
            coords1 = [corner + coord
                       for coord in self.rectangles[idx].coords]
            rects.append(polygon2d(coords1, fill=False))
        plt = sum(rects)

        if points is not None:
            for idx_rect, point in points:
                corner = corners[idx_rect]
                plt += point2d([corner + point], color='green', size=20)

        if segs is not None:
            for idx_rect, seg in segs:
                corner = corners[idx_rect]
                plt += line2d([corner + seg[0], corner + seg[1]],
                              color='orange')
        return plt

    def _get_corners(self, instrs):
        corners = [np.zeros(2) for _ in range(len(self.rectangles))]
        for idx0, dirn, idx1 in instrs:
            r0, r1 = self.rectangles[idx0], self.rectangles[idx1]
            if dirn == 'r':
                translation = np.array([r0.base, 0], dtype=type(self.field))
            elif dirn == 'l':
                translation = np.array([-r1.base, 0], dtype=type(self.field))
            elif dirn == 'u':
                translation = np.array([0, r0.height], dtype=type(self.field))
            elif dirn == 'd':
                translation = np.array([0, -r1.height], dtype=type(self.field))
            else:
                assert False
            corners[idx1] = corners[idx0] + translation
        return corners

    def decompose(self, m):
        zerovec = np.zeros(2, dtype=type(self.field))
        return [x for idx in range(len(self.rectangles))
                for x in self.flow(idx, zerovec, m)]

    def flow(self, idx_init, pt_init, m):
        data = []
        pt0, idx = np.copy(pt_init), idx_init
        while True:
            rect = self.rectangles[idx]
            pt1 = TV.step_flow(rect, pt0, m)
            data.append((idx, (pt0, pt1)))
            x1, y1 = pt1
            if idx == idx_init and TV._is_between_strict(pt0, pt_init, pt1):
                data.pop()
                data.append((idx, (pt0, pt_init)))
                return data
            elif (x1, y1) == (rect.base, rect.height):
                if self.top[self.right[idx]] != self.right[self.top[idx]]:
                    return data
                pt0 = np.zeros(2, dtype=type(self.field))
                idx = self.top[self.right[idx]]
            elif x1 == rect.base:
                pt0 = np.array([0, y1], dtype=type(self.field))
                idx = self.right[idx]
            elif y1 == rect.height:
                pt0 = np.array([x1, 0], dtype=type(self.field))
                idx = self.top[idx]
            else:
                assert False

            if idx == idx_init and (pt0 == pt_init).all():
                return data

    @staticmethod
    def step_flow(rect, pt, m):
        assert m > 0
        x, y = pt
        m_corner = (rect.height - y)/(rect.base - x)
        if 0 < m < m_corner:
            # intersects right side
            return np.array([rect.base, y + m * (rect.base - x)])
        elif m == m_corner:
            # intersects corner
            return np.array([rect.base, rect.height])
        else:
            # intersects top
            return np.array([x + (rect.height - y)/m, rect.height])

    @staticmethod
    def _is_between_strict(p0, p1, p2):
        d01 = np.linalg.norm(p1 - p0)
        d12 = np.linalg.norm(p2 - p1)
        d02 = np.linalg.norm(p2 - p0)

        return all([d01 != 0, d12 != 0, d01 + d12 == d02])
