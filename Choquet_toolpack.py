'''
Created on 13 Jan 2011

@author: gb
'''
from choquet_base import Choquet,Mobius,bitCount,Ch_gradient,MobiusB,Choquet_perm
from numpy import *
from itertools import *
from collections import defaultdict
import operator
from cap_dnf import cap_dnf
from math import factorial
import random as rnd
import cvxopt as cvx
from cvxopt import solvers
from scipy.optimize import fmin_slsqp
from time import time
from openopt import NSP
from scipy import *
from capacity_parameters import *
#from fss_test import v, zeros, dot, shape
#import networkx as nx
#import math as math
#from openopt import NSP
#cvx.solvers.options['LPX_K_MSGLEV'] = 0
### Constants
dim = 4
Fx = [lambda x: 1-exp(-3*x), lambda x: 1-exp(-3*x), lambda x: 1-exp(-3*x), lambda x: 1-exp(-3*x),lambda x: 1-exp(-3*x),lambda x: 1-exp(-x)]
#Fx = [lambda x: -0.029335929406647*pow(x,6) + 0.207346406139707*pow(x,5) - 0.579333438490725*pow(x,4) + 0.821527980327791*pow(x,3) - 0.698397454788673*pow(x,2)  + 0.687228652793793*x, lambda x: -0.027728205126421*pow(x,6) + 0.197460701672059*pow(x,5) - 0.567139344258168*pow(x,4) + 0.888084577320966*pow(x,3) - 1.002827533837975*pow(x,2)  + 1.152014732068982*x,lambda x: 1-exp(-3*x), lambda x: 1-exp(-3*x),lambda x: 1-exp(-3*x),lambda x: 1-exp(-3*x)]
#dFx = [lambda x: -0.029335929406647*6*pow(x,5) + 0.207346406139707*5*pow(x,4) - 0.579333438490725*4*pow(x,3) + 0.821527980327791*3*pow(x,2) - 0.698397454788673*2*x  + 0.687228652793793, lambda x: -0.027728205126421*6*pow(x,5) + 0.197460701672059*5*pow(x,4) - 0.567139344258168*4*pow(x,3) + 0.888084577320966*3*pow(x,2) - 1.002827533837975*2*x  + 1.152014732068982, lambda x: 3*exp(-3*x), lambda x: 3*exp(-3*x),lambda x: 3*exp(-3*x),lambda x: 3*exp(-3*x)]
#######


#def MobiusB(x):
#    """
#    "Anti-Mobius" for a non-linear case
#    """    
#    f = zeros(2**dim)
#    for v in range(1,2**dim):     # transform inverse to mobius - C = v(A) sum_{B > A} (-1)^|B\A| min_{i < B} f_i  (i.e all SUPersets)
#        for j in range(v,2**dim+1):
#            if v & j == v:
#                f_set = [p for p in range(dim) if j & (1 << p)]   # moves a bit along and checks if it is present in the integer
##                f[v] += power(-1,bitCount(j) - bitCount(v))*min([Fx[p](x[p]) for p in f_set])
#                f[v] += power(-1,bitCount(j - v))*min([Fx[p](x[p]) for p in f_set])
#                if abs(f[v]) < 0.0000001:
#                    f[v] = 0.0
#    return cvx.matrix(f)  #CVX
##    return array(f)  # Pure Mosek

def Shapley_vals(cap):
    """
    Calculates Shapley values for a given capacity
    """
    m = len(cap)
    for el in range(int(log2(m))):
        el=pow(2,el)
        mul = zeros(m)
        for i in range(m):        
            if i & el != el:
                coeff = float(factorial(bitCount(m-1) - bitCount(i) - 1))*factorial(bitCount(i))/factorial(bitCount(m-1))
                mul[i | el] = coeff
                mul[i] = -coeff
        print log2(el)+1, dot(mul,cap)
    return 0  

def MobiusC(x):
    """ 
    "Anti-Mobius" for a linear case
    """
    f = zeros(2**dim)
    for v in range(1,2**dim):     # transform inverse to mobius - C = v(A) sum_{B > A} (-1)^|B\A| min_{i < B} f_i  (i.e all SUPersets)
        for j in range(v,2**dim+1):
            if v & j == v:
                f_set = [p for p in range(dim) if j & (1 << p)]   # moves a bit along and checks if it is present in the integer
                f[v] += power(-1,bitCount(j) - bitCount(v))*min([x[p] for p in f_set])
    return cvx.matrix(f)

def find_centre():
    """
    Finds the point f_1 = ... = f_n by optimizing a capacity having v(A)=0,forall A \neq N
    """
    capacity = zeros(2**dim)
    capacity[-1] = 1
    centre = max_choquet(capacity)
    val_centre = Choquet(centre,capacity)
    return centre, val_centre

def cap_dnf_t(cap,S,cmatr = zeros((dim,dim)),result = []):
    """
    Finds the point f_1 = ... = f_n by optimizing a capacity having v(A)=0,forall A \neq N
    """
    print "IN CAPDNF"
    t0 = time()
    result = cap_dnf(cap,S,cmatr = zeros((dim,dim),dtype=int),result = [])
    print "TIME CAP_DNF", time()-t0
    print shape(result)
    return result

#Budg = 1

def max_choquet(capacity,Budg = 1):
    """
    Calculates the integral maximum
    """
    x0 = zeros(dim+1)
    def f_eqcons(x):
        Aeq = ones(len(x)-1)
        beq = Budg
        return array([dot(Aeq,x[1:])-beq])

    def fprime_eqcons(x):
        return append(0,ones(len(x)-1))

    def f_ineqcons(x,capacity):
        A = append(Choquet(x[1:],capacity)-x[0], x)
        A = append(A,Budg-x)
        return A
        

    def fprime_ineqcons(x,capacity):
        A = append(-1,Ch_gradient(x[1:], capacity))
        B = diag(ones(len(x)))        
        C = -diag(ones(len(x)))
        D = vstack((A,B,C))
        return  D

    def fprime(x):
        return append(-1,zeros(len(x)-1))

    def objfunc(x):
        return -x[0]
    
    f_ineqcons_w = lambda x: f_ineqcons(x,capacity)
    fprime_ineqcons_w = lambda x: fprime_ineqcons(x,capacity)
#    t0 = time()
#    sol = fmin_slsqp(objfunc,x0, fprime=fprime, f_eqcons=f_eqcons, f_ieqcons=f_ineqcons_w, fprime_eqcons = fprime_eqcons,iprint=2,full_output=1)
    sol = fmin_slsqp(objfunc,x0, fprime=fprime, f_eqcons=f_eqcons, f_ieqcons=f_ineqcons_w, fprime_eqcons = fprime_eqcons,fprime_ieqcons = fprime_ineqcons_w, iprint=0,full_output=1,acc=1e-11,iter=300)
    return sol[0][1:]
    
def solve_mmax(Umm,Budg = 1):
    """
    Solves mmax regret problem. Input: array of pairs (maximizer, capacity)
    """
    x0 = zeros(dim+1)
#    x0 = [0,0.1,0.2,0.3,0.4]
    
    def f_eqcons(x):
        Aeq = ones(len(x)-1)
        beq = Budg
        return array([dot(Aeq,x[1:])-beq])
    
    def fprime_eqcons(x):
        return append(0,ones(len(x)-1))
    
    def f_ineqcons(x,Umm):
        A = [x[0] - (Choquet(array(p),v) - Choquet(x[1:],v)) for p,v in Umm]
        A.extend(x)
#        A.extend(x[1:])
        return array(A)
    
    def f_ineqcons_prime(x,Umm):
        A = array([append(1,Ch_gradient(x[1:], v)) for p,v in Umm])
        B = diag(ones(len(x)))        
        D = vstack((A,B))
        return  D
    
    def fprime(x):
        return append(1,zeros(len(x)-1))

    def objfunc(x):
        return x[0]
    
    f_ineqcons_w = lambda x: f_ineqcons(x,Umm)
    fprime_ineqcons_w = lambda x: f_ineqcons_prime(x,Umm)
#    print capacity
    sol = fmin_slsqp(objfunc,x0, fprime=fprime, f_eqcons=f_eqcons, f_ieqcons=f_ineqcons_w, fprime_eqcons = fprime_eqcons, fprime_ieqcons = fprime_ineqcons_w, iprint=2,full_output=1,acc=1e-11,iter=900)
#    sol = fmin_slsqp(objfunc,x0, fprime=fprime, f_eqcons=f_eqcons, f_ieqcons=f_ineqcons_w, fprime_eqcons = fprime_eqcons, iprint=2,full_output=1,acc=1e-11,iter=900)
    print sol
    return sol[0][1:],sol[1]

########## Tree

def node_val(G,n,f_vals):
    """
    G - graph, n - the root node. Calculates its value based on values of all children
    """
    G_pred = G.node[n]['psorted']
    if not G_pred:
        f_weights = defaultdict(lambda: 0)
        f_weights[n] = 1
        return G.node[n]['func'](f_vals[n]), f_weights 
    else:
        X = [node_val(G,p,f_vals) for p in G_pred]
        n_weights = Choquet_perm(array([p[0] for p in X]),G.node[n]['cap'])
        f_weights = dot(array([[p[1][q] for p in X] for q in f_vals.keys()]),n_weights)
#        print [G.node[i]['func'](f_vals[i]) for i in f_vals]
        n_value = dot([G.node[i]['func'](f_vals[i]) for i in f_vals],f_weights)
        f_weights = defaultdict(lambda: 0,zip(f_vals.keys(),f_weights))
#        print n, n_value
#        print f_weights
        return n_value,f_weights 


def max_tree(G,n,fsorted,Budg = 1, x0 = array([])):
    """
    G - graph, n - the root node. Maximizes the function value at node n, distributing Budg between children
    """
    if not x0.shape[0]:
        x0 = zeros(len(fsorted)+1)

    def f_eqcons(x):
        Aeq = ones(len(x)-1)
        beq = Budg
        return array([dot(Aeq,x[1:])-beq])

    def fprime_eqcons(x):
        return append(0,ones(len(x)-1))

    def f_ineqcons(x,G,n,fsorted):
        n_value = node_val(G,n,dict(zip(fsorted,x[1:])))[0]
        A = append(n_value-x[0], x)
        A = append(A,Budg-x)
#        A = append(A,1-x)
        return A
        

    def fprime_ineqcons(x,G,n,fsorted):
        f_val = dict(zip(fsorted,x[1:]))
        weights = node_val(G,n,f_val)[1]
        grad = [G.node[p]['gradient'](f_val[p])*weights[p] for p in fsorted]
        A = append(-1,grad)
        B = diag(ones(len(x)))        
        C = -diag(ones(len(x)))
        D = vstack((A,B,C))
        return  D

    def fprime(x):
        return append(-1,zeros(len(x)-1))

    def objfunc(x):
        return -x[0]
    
    f_ineqcons_w = lambda x: f_ineqcons(x,G,n,fsorted)
    fprime_ineqcons_w = lambda x: fprime_ineqcons(x,G,n,fsorted)
#    t0 = time()
#    sol = fmin_slsqp(objfunc,x0, fprime=fprime, f_eqcons=f_eqcons, f_ieqcons=f_ineqcons_w, fprime_eqcons = fprime_eqcons,iprint=2,full_output=1,epsilon=4e-08)
    sol = fmin_slsqp(objfunc,x0, fprime=fprime, f_eqcons=f_eqcons, f_ieqcons=f_ineqcons_w, fprime_eqcons = fprime_eqcons,fprime_ieqcons = fprime_ineqcons_w, iprint=1,full_output=1,acc=1e-13,iter=3000)
#    a = time()-t0
#    print a, sol[1]
#    print sol
#    if a > 2:
#        print capacity
    return sol[0][1:],sol[0][0]

###SUBGRAD

def simpl_project(x):
    """
    Simplex projection algorithm 
    Shalev-Shwartz, S., & Singer, Y. (2006). Efficient learning of label ranking by soft projections onto polyhedra
    Duchi, J. and Shalev-Shwartz, S. and Singer, Y. and Chandra, T. Efficient projections onto the l1-ball for learning in high dimensions
    """
    m = sorted(x,reverse=1)
    r = max([j for j in range(1,len(m)+1) if (m[j-1] - 1./j*(sum(m[0:j-1]) - 1))  > 0])
    th = 1./r*(sum(m[0:r])-1)
    return array([max(p-th,0) for p in x])

def max_choquet_sbgr(G,n,fsorted):
    """
    subgradient projection maximization. Expects graph, node, and list of grpah functions as input
    """
    x =zeros(len(fsorted))
    iter = 10000
    fbest = 0
    xbest = 0
    for i in range(iter):
        if (x < 0).any():
            e = zeros(len(fsorted))
            e[nonzero(x<0)[0]] = -1
        elif (x > 1).any():
            e = zeros(len(fsorted))
            e[nonzero(x>1)[0]] = 1
        else:
            f_val = dict(zip(fsorted,x))
            weights = node_val(G,n,f_val)[1]
            grad = [G.node[p]['gradient'](f_val[p])*weights[p] for p in fsorted]
            e = -array(grad)
        x = simpl_project(x - pow((i+1),-0.9)*e)
        f = node_val(G,n,dict(zip(fsorted,x)))[0]
        if f > fbest:
            fbest = f
            xbest = x
    return fbest, xbest

def max_dual(cap_array,Budg = 1):
    """
    Theoretically calculates min(z^r)->max(lambda) i.e. the dual problem solution, but in fact does not work
    """
    x0 = array([1./shape(cap_array)[0] for i in range(shape(cap_array)[0])])
    print x0
    def f_eqcons(x):
#        Aeq = ones(len(x))
        beq = 1
        print "ARRR", array([sum(x)-beq])
        return array([sum(x)-beq])

    def fprime_eqcons(x):
        return -ones(len(x))

    def f_ineqcons(x):
        return x
        
    def fprime_ineqcons(x):
        #A = append(-1,Ch_gradient(x[1:], capacity))
        B = diag(ones(len(x)))        
        #C = -diag(ones(len(x)))
        #D = vstack((A,B,C))
        return  B

    def fprime(x,cap_array,Budg):
        cap_int = dot(cap_array.T,x)
        x_int = array(max_choquet(cap_int,Budg))
        ch_diffs = array([Choquet(x_int,v) - Choquet(array(max_choquet(v,Budg)),v) for v in cap_array])
        print "GRAD"
        return ch_diffs

    def objfunc(x,cap_array,Budg):
        cap_int = dot(cap_array.T,x)
        x_int = array(max_choquet(cap_int,Budg))
        ch_diffs = array([Choquet(x_int,v) - Choquet(array(max_choquet(v,Budg)),v) for v in cap_array])
        print dot(x, ch_diffs) 
        return dot(x, ch_diffs)
    
    objfunc_w = lambda x: objfunc(x,cap_array,Budg)
    fprime_w = lambda x: fprime(x,cap_array,Budg)
#    sol = fmin_slsqp(objfunc_w,x0, fprime=fprime_w, f_eqcons=f_eqcons, f_ieqcons=f_ineqcons, iprint=2,full_output=1,iter=300)
    sol = fmin_slsqp(objfunc_w,x0, fprime=fprime_w, f_eqcons=f_eqcons, f_ieqcons=f_ineqcons, fprime_eqcons = fprime_eqcons,fprime_ieqcons = fprime_ineqcons, iprint=2,full_output=1,acc=1e-11,iter=300)
    print sol
    return sol[0][1:]


def max_choquet_glob(capacity):
    """
    Calculates the global maximum of Choq. w.r.t any capacity by decomposing into convex parts
    """
    mcap = Mobius(capacity)
    sols = [[max_choquet(p),p] for p in cap_dnf_t(mcap,nonzero(mcap<0)[0][0],cmatr = zeros((dim,dim),dtype=int),result = [])]
    s1 = [[Choquet(array(p[0]),p[1]),p[0]] for p in sols]
    return max(s1, key=operator.itemgetter(0))[1]
    
def solve_mmax_wval(Umm,Budg = 1):
    """
    In this version Umm is suppposed to hold a VALUE in [0]s, not the point
    """
    x0 = zeros(dim+1)
#    x0 = [0,0.1,0.2,0.3,0.4]
    
    def f_eqcons(x):
        Aeq = ones(len(x)-1)
        beq = Budg
        return array([dot(Aeq,x[1:])-beq])
    
    def fprime_eqcons(x):
        return append(0,ones(len(x)-1))
    
    def f_ineqcons(x,Umm):
        A = [x[0] - (p - Choquet(x[1:],v)) for p,v in Umm]
        A.extend(x)
#        A.extend(x[1:])
        return array(A)
    
    def f_ineqcons_prime(x,Umm):
        A = array([append(1,Ch_gradient(x[1:], v)) for p,v in Umm])
        B = diag(ones(len(x)))        
        D = vstack((A,B))
        return  D
    
    def fprime(x):
        return append(1,zeros(len(x)-1))

    def objfunc(x):
        return x[0]
    
    f_ineqcons_w = lambda x: f_ineqcons(x,Umm)
    fprime_ineqcons_w = lambda x: f_ineqcons_prime(x,Umm)
#    print capacity
    sol = fmin_slsqp(objfunc,x0, fprime=fprime, f_eqcons=f_eqcons, f_ieqcons=f_ineqcons_w, fprime_eqcons = fprime_eqcons, fprime_ieqcons = fprime_ineqcons_w, iprint=2,full_output=1,acc=1e-11,iter=900)
#    sol = fmin_slsqp(objfunc,x0, fprime=fprime, f_eqcons=f_eqcons, f_ieqcons=f_ineqcons_w, fprime_eqcons = fprime_eqcons, iprint=2,full_output=1,acc=1e-11,iter=900)
    print sol
    return sol[0][1:],sol[1]


def solve_dual(cap_array,gm):
    """
    Input: array of tuples (capacity (e.g. vertices of U) or its nec.measure component, global maximum)
    Finds the robust capacity, i.e. the solution of dual problem
    Capacities assumed to be 2-MONOTONE
    Output: v^r
    """
    Budg = 1  
    x0 = array([1./shape(cap_array)[0] for i in range(shape(cap_array)[0])])
#    print x0
    Aeq = ones(shape(cap_array)[0])
    beq = 1
    A = diag(-ones(shape(cap_array)[0]))
    b = zeros(shape(cap_array)[0])
    vmax=gm
    # vmax = array([chq.Choquet(array(chq.max_choquet(v,Budg)),v) for v in cap_array])
    
    def fprime(x,cap_array,Budg):
        cap_int = dot(cap_array.T,x)
        x_int = array(max_choquet(cap_int,Budg))
        ch_xr = array([Choquet(x_int,v) for v in cap_array])
        ch_diffs = ch_xr - vmax
        return ch_diffs

    def objfunc(x,cap_array,Budg):
        cap_int = dot(cap_array.T,x)
        x_int = array(max_choquet(cap_int,Budg))
        ch_xr = array([Choquet(x_int,v) for v in cap_array])
        ch_diffs = ch_xr - vmax
#        print  dot(x, ch_diffs)
        return dot(x, ch_diffs)
    
    objfunc_w = lambda x: objfunc(x,cap_array,Budg)
    fprime_w = lambda x: fprime(x,cap_array,Budg)
    
    p = NSP(objfunc_w, x0, df=fprime_w,  A=A,  b=b,  Aeq=Aeq,  beq=beq, iprint = 25, maxIter = 300, maxFunEvals = 1e7, diffint=1e-10, xtol=1e-9,ftol=1e-9,gtol=1e-5)
    r = p.solve('algencan')
    print r.xf, r.ff
    cap_mix = dot(cap_array.T, array(r.xf))
    return r.ff,r.xf
