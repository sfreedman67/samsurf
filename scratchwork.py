import random
import more_itertools
from more_itertools import first_true

A = [1, 2, 3, 4, 20, 21, 22, 23, 5, 6, 7, 8, 9]

# Want [1, 2, 3, 4] + [NEW POINTS] + [5, 6, 7, 8, 9]
# edge_intersected_first: ai = T but ai+1 = F
# since the true region is on your left and polygon oriented CCW
# edge_intersected_second : ai = F but ai+1 = T
# So, want to return idx(4) = 3 and idx(23) = 7 in this example
# since then we can write A[:3 + 1] for the first batch of trues, and A[7 + 1:] for the second
# THe +1 is awkard, What if we return start of the Falses, start of the second Trues


edges = list(zip(A, A[1:]+A[:1]))
print(list(edges))
print(first_true(range(len(edges)), pred=lambda i: A[i] < 10 and A[i + 1] > 10))

