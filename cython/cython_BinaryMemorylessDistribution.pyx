# from math import log2

cdef extern from "math.h":
        double log2(double x) nogil

# functions for degrade/upgrade/merge
cpdef float eta(float p):
    if p == 0:
        return 0
    else:
        return -p * log2(p)

cpdef float hxgiveny(data):
    cdef float d0 = data[0]
    cdef float d1 = data[1]
    
    cdef float py = d0 + d1
    return py * ( eta(d0/py) + eta(d1/py) )
