from triangulation import Triangle, Triangulation
from sage.all import *
import flatsurf as fs
import typing


# can I really keep flatsurf methods contained here?
def dlny_triang(X: fs.TranslationSurface) -> Triangulation:
    DT = X.delaunay_triangulation()

    tris = []

    for i in range(DT.num_polygons()):
        tri = Triangle([vector(edge) for edge in DT.polygon(i).edges()])
        tris.append(tri)

    gluings = {edge[0]: edge[1] for edge in DT.edge_iterator(gluings=True)}

    return Triangulation(tris, gluings)


if __name__ == "__main__":
    # eventually, user will select their surface
    X = fs.translation_surfaces.regular_octagon()
    DT = dlny_triang(X)

    g = HyperbolicPlane().UHP().get_geodesic(2,3)
    # TODO: how to show here?
    print(g.plot(axes=True))


    # Want to make a test case for the edge ineqs and sage arithmetic
    print(DT.triangle(0))
