# N = number of triangles , numbered 0 \le i < N
# j = 0,1,2 denotes side of triangle
# Must have 2E = 3T

#An edge is of the form (i, j), where i is a triangle label, j an edge label 
#A triangulation is a *DICT* of gluings e1 -> e2, where e1, e2 are edges of two different triangles.
#This has the annoying feature that if e1:e2 is an entry, then e2:e1 also needs to be a valid entry
#Why use set for a triangulation? I want to be able to tell quickly if two edges are glued
#Validity = every edge has a unique glued partner on a different triangle
#Ex of a triangulation: {((1, 0), (2, 1)),...] means triangle 1, side 0 glues to triangle 2, side 1

#A combinatorial equivalence consists of a bijection on the triangles that respects edge gluings
#When we biject triangles, we need a "shift" to tell how we glued

#Brute-force algorithm to find all equivalences  
#Step 1: generate all possible "bijections" (perm + shift) on triangles
#Step 2: see if each gluing respects the edge identifications


import itertools as it
import unittest
from sage.all import *
import flatsurf
import sage.graphs.generic_graph as sg
from collections import deque

def set_edge_gluings(t):
	return set(tuple(sorted(gluing)) for gluing in t.edge_iterator(gluings=True))

def gen_poss_homeos(N):
	return set(it.product(it.permutations(range(N)),it.product(range(3),repeat=N)))

def triang_to_graph(t):
	G = Graph(multiedges=True)
	
	for edge_pair in set_edge_gluings(t):
			edge1 = edge_pair[0]
			edge2 = edge_pair[1]
			tri1 = edge1[0]
			tri2 = edge2[0]
			G.add_edge(tri1, tri2)

	return G

#Returns where an edge maps under a comb equiv
def image_edge(equiv, edge):
	perm = equiv[0]
	shifts = equiv[1]

	tri_lab = edge[0]
	edge_lab = edge[1]

	return (perm[tri_lab], (edge_lab + shifts[tri_lab]) % 3)

def homeo_respects_gluings(equiv, t1, t2):
	#for all gluings (e1,e2) in t1, want to see if F(e1) and F(e2) are glued in t2
	for e1 in t1.edge_iterator():
		e2 = t1.opposite_edge(e1)
		f1 = image_edge(equiv, e1)
		f2 = image_edge(equiv, e2)
		if t2.opposite_edge(f1) != f2:
			return False

	return True

def rotate_edge(e, k):
	return (e[0], (e[1] + k) % 3)

def neighbor_inferences(e,f):
	return {rotate_edge(e, k): rotate_edge(f, k) for k in range(3)}

def edges_opp_tri(t, n):
	return [t.opposite_edge(n,k) for k in range(3)]

#extends an assignment (0,0) -> f to a mapping. It *might not* be well-def or CE
def generate_homeo_from_edge(t1, t2, f):
	num_tris = t1.num_polygons()

	#tree search through the triangles. We enter a triangle at its base side, exit through the other two 
	
	edges_to_visit = deque(edges_opp_tri(t1, 0))

	tris_seen = {0}

	partial_mapping = neighbor_inferences((0,0), f)

	#Check against BFS pseudocode
	while(len(tris_seen) < num_tris):
		curr_edge = edges_to_visit.pop()
		curr_tri = curr_edge[0]
		#we assume that we got to curr_edge from its opposite edge
		if curr_tri not in tris_seen: 
			opp_edge = t1.opposite_edge(curr_edge)
			image_opp_edge = partial_mapping[opp_edge]
			image_curr_edge = t2.opposite_edge(image_opp_edge) 
			partial_mapping.update( neighbor_inferences(curr_edge, image_curr_edge))
			tris_seen.add(curr_tri)

		#we now exit out the two neighbors		
		edges_to_visit.extendleft(edges_opp_tri(t1, curr_tri))

	#now, we want to turn out partial mapping into a homeo
	perm = []
	shifts = [] 
	
	for n in range(num_tris):
		#for each n, find an edge on triangle n
		edges_n = [edge for edge in partial_mapping if edge[0] == n]
		assert(edges_n)
		edge_n = edges_n[0]
		image_edge_n = partial_mapping[edge_n]

		perm.append(image_edge_n[0])
		shifts.append((image_edge_n[1] - edge_n[1]) % 3)

	return (tuple(perm), tuple(shifts))


#Returns a list of all combinatorial equivalences 
def gen_comb_equivs(t1, t2):
	# brute-force method was: set(filter(lambda x: homeo_respects_gluings(x,t1,t2), gen_poss_homeos(N)))
	assert(t1.num_polygons() == t2.num_polygons())
	N = t1.num_polygons()

	poss_homeos = {generate_homeo_from_edge(t1, t2, edge) for edge in it.product(range(N), range(3))}
	return set(filter(lambda x: homeo_respects_gluings(x, t1, t2), poss_homeos))


class MyTests(unittest.TestCase):
	def test_generate_homeo_from_edge(self):
		t = flatsurf.translation_surfaces.square_torus().delaunay_triangulation()

		self.assertEqual(generate_homeo_from_edge(t,t,(0,0)), ((0,1), (0,0)))
		self.assertEqual(generate_homeo_from_edge(t,t,(1,1)), ((1,0),(1,1)))

	def test_gen_poss_homeos(self):
		self.assertEqual(len(gen_poss_homeos(2)), 18)
		self.assertEqual(len(gen_poss_homeos(6)), 524880)

	def test_torus_self_equivs(self):
		t = flatsurf.translation_surfaces.square_torus()
		DT_sq_torus = t.delaunay_triangulation()
		equivs_two_tris = set(it.product([(1,0),(0,1)], [(i,i) for i in range(3)]))
		self.assertEqual(gen_comb_equivs(DT_sq_torus,DT_sq_torus), equivs_two_tris)

	def test_veech_octagon_no_comb_autos(self):
		v_oct = flatsurf.translation_surfaces.regular_octagon()
		DT_v_oct = v_oct.delaunay_triangulation()
		self.assertEqual(gen_comb_equivs(DT_v_oct,DT_v_oct), {((0,1,2,3,4,5),(0,0,0,0,0,0))})

unittest.main(verbosity=2)