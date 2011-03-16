'''
Created on 13 Jan 2011

@author: gb
'''
from choquet_base import Choquet,Mobius,bitCount,Ch_gradient,MobiusB,Choquet_perm
from numpy import *
from itertools import *
from collections import defaultdict
import networkx as nx
from cap_dnf import cap_dnf
#from math import power
#import math as math
import random as rnd
import cvxopt as cvx
from cvxopt import solvers
from scipy.optimize import fmin_slsqp
#from openopt import NSP
from time import time
#cvx.solvers.options['LPX_K_MSGLEV'] = 0
### Constants
dim = 3
Fx = [lambda x: 1-exp(-3*x), lambda x: 1-exp(-3*x), lambda x: 1-exp(-3*x), lambda x: 1-exp(-3*x),lambda x: 1-exp(-3*x),lambda x: 1-exp(-3*x)]
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
    Finds the centre point
    """
    capacity = zeros(2**dim)
    capacity[-1] = 1
    centre = max_choquet(capacity)
    val_centre = Choquet(centre,capacity)
    return centre, val_centre

def cap_dnf_t(cap,S,cmatr = zeros((dim,dim)),result = []):
    print "IN CAPDNF"
    t0 = time()
    result = cap_dnf(cap,S,cmatr = zeros((dim,dim),dtype=int),result = [])
    print "TIME CAP_DNF", time()-t0
    print shape(result)
    return result
"""
def max_choquet_00(capacity):
#    print capacity
    x0 = [0.2, 0.3, 0.5]
    A = -eye(dim)
    b = zeros(dim)
    Aeq = ones(dim)
    beq = 1
    cwrap_vx = lambda x: Choquet(x,capacity)
    x_opt = NSP(cwrap_vx, x0, A=A, b=b, Aeq=Aeq, beq=beq)
    r = x_opt.maximize('ralg')
    return r.xf
"""
#Budg = 1

def max_choquet(capacity,Budg = 1.5):
#    x0 = [0,0.2,0.3,0.5]
    x0 = zeros(dim+1)
#    x0 = [0,1./3,1./3,1./3]
#    x0 = [1,0.,0,1]
#    x0 = [1,0.25,0.25,0.25,0.25]
#    x0 = [1,0.2,0.2,0.2,0.2,0.2]
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
#    a = time()-t0
#    print a, sol[1]
#    print sol
#    if a > 2:
#        print capacity
    return sol[0][1:]
    
def solve_mmax(Umm,Budg = 1.5):
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
#        A = append(A,Budg-x)
        A = append(A,1-x)
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
    sol = fmin_slsqp(objfunc,x0, fprime=fprime, f_eqcons=f_eqcons, f_ieqcons=f_ineqcons_w, fprime_eqcons = fprime_eqcons,fprime_ieqcons = fprime_ineqcons_w, iprint=2,full_output=1,acc=1e-13,iter=3000)
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