# from math import log2

cdef extern from "math.h":
        double log2(double x) nogil

# functions for degrade/upgrade/merge
cpdef float eta(float p):
    if p == 0:
        return 0
    else:
        return -p * log2(p)

cpdef double hxgiveny(list data):
    cdef double d0 = data[0]
    cdef double d1 = data[1]
    
    cdef double py = d0 + d1
    return py * ( eta(d0/py) + eta(d1/py) )
