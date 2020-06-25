""" Given an arbitrary q-triangulation of X, return a Delaunay triangulation of X.
usage::
	>>> import edge_flip, triangulation
	>>> T = foo_triangulation
	>>> DT = edge_flip(T)
:param T: A Triangulation (i.e. triangles and edge gluings) of an initial flat surface.
:rtype: A Triangulation, that is Delaunay  
"""
