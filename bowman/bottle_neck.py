import sage.all
from sage.all import *


def _is_negative(A, B, C):
    if B == 0 or C == 0:
            return bool(A < 0)
    K = (-A) / B
    if B > 0:
        return bool(K >= 0) and bool((C - K**2) < 0)
    return K <= 0 or C - K**2 > 0

def _is_zero(A, B, C):
    if B == 0 or C == 0:
        return A == 0
    
    K = (-A) / B
    return bool(K >= 0) and bool(C - K**2 == 0)


def _bottle_neck(A, B, C, only_zero):
    if only_zero:
        return _is_zero(A, B, C)
    return not _is_negative(A, B, C)     
        

                
