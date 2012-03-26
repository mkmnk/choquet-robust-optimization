cimport cython
import numpy as np
cimport numpy as np
from itertools import permutations
from random import shuffle
#from time import time
cdef extern from "math.h":
  double pow(double,double)
  double exp(double)

#cdef extern from "time.h":
#   ctypedef long time_t
#   time_t time(time_t*)


Fx = [lambda x: 1-exp(-3*x), lambda x: 1-exp(-3*x), lambda x: 1-exp(-3*x), lambda x: 1-exp(-3*x),lambda x: 1-exp(-3*x),lambda x: 1-exp(-x)]
dFx = [lambda x: 3*exp(-3*x), lambda x: 3*exp(-3*x), lambda x: 3*exp(-3*x), lambda x: 3*exp(-3*x),lambda x: 3*exp(-3*x),lambda x: exp(-x)]
#Fx = [lambda x: -0.029335929406647*pow(x,6) + 0.207346406139707*pow(x,5) - 0.579333438490725*pow(x,4) + 0.821527980327791*pow(x,3) - 0.698397454788673*pow(x,2)  + 0.687228652793793*x, lambda x: -0.027728205126421*pow(x,6) + 0.197460701672059*pow(x,5) - 0.567139344258168*pow(x,4) + 0.888084577320966*pow(x,3) - 1.002827533837975*pow(x,2)  + 1.152014732068982*x,lambda x: 1-exp(-3*x), lambda x: 1-exp(-3*x),lambda x: 1-exp(-3*x),lambda x: 1-exp(-3*x)]
#dFx = [lambda x: -0.029335929406647*6*pow(x,5) + 0.207346406139707*5*pow(x,4) - 0.579333438490725*4*pow(x,3) + 0.821527980327791*3*pow(x,2) - 0.698397454788673*2*x  + 0.687228652793793, lambda x: -0.027728205126421*6*pow(x,5) + 0.197460701672059*5*pow(x,4) - 0.567139344258168*4*pow(x,3) + 0.888084577320966*3*pow(x,2) - 1.002827533837975*2*x  + 1.152014732068982, lambda x: 3*exp(-3*x), lambda x: 3*exp(-3*x),lambda x: 3*exp(-3*x),lambda x: 3*exp(-3*x)]

@cython.boundscheck(False)
@cython.wraparound(False)

cdef inline int i_min(int a, int b): return a if a <= b else b
cdef inline double d_max(double a, double b): return a if a >= b else b

cdef np.ndarray[np.int_t, ndim=2] WFL(np.ndarray[np.int_t, ndim=2] A):
    """ 
    Warshall-Floyd 
    """
#    cdef float t0 = time()
    cdef int dim = 0 
    cdef int k = 0
    cdef int i = 0
    cdef int j = 0
    cdef int infty = 9999
    np.putmask(A,A == 0,infty)
    dim = A.shape[0]
    for k in range(dim):
        for i in range(dim):
            for j in range(dim):
                A[i,j] = i_min(A[i,j], (A[i,k] + A[k,j]))
    np.putmask(A,A == infty,0)
#    print "WFL TIME", t0-time()
    return A

cdef np.ndarray[np.float64_t, ndim=1] zeta(np.ndarray[np.float64_t, ndim=1] cap):
    """
    Inverse to mobius
    """
    cdef np.ndarray[np.float64_t, ndim=1] nu = np.zeros(len(cap))
    cdef int i = 0
    cdef int j = 0
    for i in range(1,len(cap)):
        for j in range(0,i+1):
            if i & j == j:
                nu[i] += cap[j]
    return nu

cdef np.ndarray[np.float64_t, ndim=1] to_bel(np.ndarray[np.float64_t, ndim=1] cap,np.ndarray[np.int_t, ndim=2] A):
    """
    Generates a belief function (bT) from set of constraints (R subset)
    """
#    cdef float t0 = time()
    cdef np.ndarray[np.float64_t, ndim=1] capnew = np.zeros(len(cap))
    cdef int i = 0
    cdef int p = 0
    cdef int j = 0
    cdef int k = 0
    for i in range(1,len(cap)):
        el_i = [p for p in range(np.binary_repr(i).__len__()) if i & (1 << p)]
        for j in el_i:
#if np.array_equal(np.intersect1d(np.where(A[j,:] > 0)[0],el_i), np.where(A[j,:] > 0)[0]):
            if set(np.where(A[j,:] > 0)[0]).intersection(set(el_i)) == set(np.where(A[j,:] > 0)[0]):
                capnew[i] = cap[i]
            else:
                capnew[i] = 0
                for k in range(0,i+1):
                    if k & i == k:
                      capnew[i] = d_max(capnew[i],capnew[k])
                break
#    print "TO_BEL TIME", t0 - time()            
    return capnew


def cap_dnf(np.ndarray[np.float64_t, ndim=1] cap, int S,np.ndarray[np.int_t, ndim=2] cmatr ,list result = []):
    """
    Breaks a non-convex capacity into convex subsets. Employs splitting algo from 
    M. Timonin, Maximization of the Choquet integral over a convex set and its application to resource allocation problems, Annals of Operations Research.
    """
    cdef int i = 0
    cdef int p = 0
    cdef int j = 0
    cdef int el1
    cdef int el2
    cdef int newin

#    cdef t0 = time()
    S_el = [p for p in range(np.binary_repr(S).__len__()) if S & (1 << p)]
#    shuffle(S_el)
#    print S_el
    for el in S_el:
        A = np.array(cmatr)
        for i in S_el:
          if i != el:
            A[el,i] = 1
        A = WFL(A)
        capnew = np.array(cap)
### Rewrite later for performance 
        for i in range(1,len(cap)):
            el_i = [p for p in range(np.binary_repr(i).__len__()) if i & (1 << p)]
            new_el = el_i[:]
            for j in el_i: 
                inf_set = set(el_i).intersection(set(np.nonzero(A[j,:])[0]))
                new_el = set(new_el).difference(set(inf_set))
            newin = sum([pow(2,p) for p in new_el])
            if newin != i:
                capnew[newin] = (capnew[newin] + cap[i])
                if (capnew[newin] > -0.0000001) and (capnew[newin] < 0.0000001):
                    capnew[newin]=0.0
                capnew[i] = 0.0
        if np.nonzero(capnew < -0.0000001)[0].any():
            capnew = cap_dnf(capnew, np.nonzero(capnew < 0)[0][0], A,result)
        else: 
            result.append(to_bel(zeta(capnew),A))
#    print "DNF ALL", t0-time()
    return  result
