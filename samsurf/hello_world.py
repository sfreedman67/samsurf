from context import samsurf
from samsurf import triangulation

X = triangulation.Triangulation.ronen_l(44)
print("This is some info about an L-shaped table")
print(X.triangles)
print(X.gluings)