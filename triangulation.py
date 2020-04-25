from sage.all import *

class Triangle:

    def __init__(self, edges):
        if edges[0] + edges[1] != -edges[2]:
            raise ValueError("sides are oriented incorrectly")
        else:
            self.edges = edges

    def edge(self, n):
        return self.edges[n]

    def __repr__(self):
        return str(self.edges)

    def __str__(self):
        return str(self.edges)


class Triangulation:

    def __init__(self, triangles, gluings):
        self.triangles = triangles
        self.gluings = gluings

    def triangle(self, n):
        return self.triangles[n]

    def edge_set(self):
    	return None

if __name__ == "__main__":
    T1 = Triangle([vector([1, 0]), vector([-1, 1]), vector([0, -1])])
    T2 = Triangle([vector([0, 1]), vector([-1, 0]), vector([1, -1])])
    gluings = {(0, 0): (1, 0), (0, 1): (1, 1), (0, 2): (1, 2),
       	       (1, 0): (0, 0), (1, 1): (0, 1), (1, 2): (0, 2)}
    Torus = Triangulation([T1, T2], gluings)
    assert(Torus.edge_set() == {(0, 0), (0, 1), (0, 2)})
