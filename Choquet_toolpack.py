'''
Created on 13 Jan 2011

@author: gb
'''
from choquet_base import Choquet,Mobius,bitCount,Ch_gradient,MobiusB
from numpy import *
from itertools import *
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
dim = 5
Fx = [lambda x: 1-exp(-3*x), lambda x: 1-exp(-3*x), lambda x: 1-exp(-3*x), lambda x: 1-exp(-3*x),lambda x: 1-exp(-3*x),lambda x: 1-exp(-3*x)]
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
#    sol = fmin_slsqp(objfunc,x0, fprime=fprime, f_eqcons=f_eqcons, f_ieqcons=f_ineqcons_w, fprime_eqcons = fprime_eqcons,iprint=2,full_output=1,epsilon=4e-08)
    sol = fmin_slsqp(objfunc,x0, fprime=fprime, f_eqcons=f_eqcons, f_ieqcons=f_ineqcons_w, fprime_eqcons = fprime_eqcons,fprime_ieqcons = fprime_ineqcons_w, iprint=0,full_output=1,acc=1e-09)
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
        A.extend(x[1:])
        return array(A)
    
    def fprime(x):
        return append(1,zeros(len(x)-1))

    def objfunc(x):
        return x[0]
    
    f_ineqcons_w = lambda x: f_ineqcons(x,Umm)
    #fprime_ineqcons_w = lambda x: fprime_ineqcons(x,cap)
#    print capacity
    sol = fmin_slsqp(objfunc,x0, fprime=fprime, f_eqcons=f_eqcons, f_ieqcons=f_ineqcons_w, fprime_eqcons = fprime_eqcons,iprint=1,full_output=1,acc=1e-09)
    return sol[0][1:],sol[1]
