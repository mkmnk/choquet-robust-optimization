'''
Created on 16 Dec 2010

@author: gb
'''

from numpy import *
from itertools import *
from math import factorial
import random as rnd
import cvxopt as cvx
from cvxopt import solvers
from time import time
import Choquet_toolpack as chq
cvx.solvers.options['show_progress'] = False
#import psyco
#psyco.full()
set_printoptions(edgeitems='Inf',linewidth=90)
#cvx.printing.options['dformat'] = '%.3f'
cvx.printing.options['width'] = -1
cvx.solvers.options['LPX_K_MSGLEV'] = 0
#### Parameters
Shapley = [(0b00001,0.2),(0b00010,0.3),(0b00100,0.1),(0b01000,0.2)]
#Shapley = [(0b00001,0.3),(0b00010,0.3),(0b00100,0.1)]
#Shapley = [(0b00001,0.5),(0b00010,0.3),(0b00100,0.1)]
#Shapley = [(0b00001,0.2),(0b00010,0.3),(0b00100,0.1)]
#Shapley = [(0b00001,0.7),(0b00010,0.1)]
#Shapley = [(0b00001,0.2),(0b00010,0.3)]
#Shapley = [(0b00100,0.8)]
#Shapley = [(0b00010,0.3)]
#Shapley = [(0b00001,0.2)]
maxiter = 10
#####
#
######MAIN

def convexity(m):
    x = []
    I = []
    J = []
    j = 0
    for i in range(1,int(log2(m))):
        combs = combinations([p for p in range(1,m) if chq.bitCount(p) == i ],2)      # change to get different k-monotonicity
        for (seta,setb) in combs:
                x.extend([1,1,-1,-1])
                I.extend([j,j,j,j])
                J.extend([seta | setb, seta & setb, seta, setb])
                j = j+1
    A = cvx.spmatrix(x,I,J)                    
    return A   

def limits(m):
    x = ones(m-2)
    I = range(m-2)
    J = range(1,m-1)   
    A = cvx.spmatrix(x,I,J,(m-2,m))
    return A

def equalities(m):
    A = cvx.spmatrix([1,1],[0,1],[0,m-1],(2,m))
    return A

def shapley(m,el):
    x = []
    I = []
    J = []
    for i in range(m):        
        if i & el != el:
            coeff = float(factorial(chq.bitCount(m-1) - chq.bitCount(i) - 1))*factorial(chq.bitCount(i))/factorial(chq.bitCount(m-1))
            x.extend([coeff,-coeff])
            I.extend([0,0])   # somehow it wouldn't take I=zeros(m)
            J.extend([i | el, i])
    A = cvx.spmatrix(x,I,J)
    return A  

def gen_inequalities(dim):
    A = cvx.sparse([-convexity(2**dim),-limits(2**dim)])
    b = cvx.matrix(0,(A.size[0],1),"d")
    return A,b

def gen_equalities(dim,Shpl):
    Alist = [equalities(2**dim)]
    Alist.extend([shapley(2**dim,p[0]) for p in Shpl])
    Aeq = cvx.sparse(Alist)
    blist = [0,1]
    blist.extend([p[1] for p in Shpl])
    beq = cvx.matrix(blist)
    return Aeq,beq

def xr_discr(xr, threshold, A, b, Aeq, beq):
#    print xr
#    print chq.MobiusB(xr).T
#    raw_input("Press Enter to continue...")
    Aeq = cvx.sparse([Aeq,chq.MobiusB(xr).T])
    beq = cvx.matrix([beq,threshold])
    zenith = zeros(size(A)[1])
    for i in range(1,size(A)[1]):
        f = cvx.matrix(zeros(size(A)[1]))
        f[i] = -1
        sol  = cvx.solvers.lp(f,A,b,Aeq,beq,'glpk')
        zenith[i] = sol['x'][i]
    return zenith

def xr_scan(xr, A, b, Aeq, beq):
    xr_max = -cvx.solvers.lp(-chq.MobiusB(xr),A,b,Aeq,beq,'glpk')['primal objective']
    xr_min = cvx.solvers.lp(chq.MobiusB(xr),A,b,Aeq,beq,'glpk')['primal objective']
    if xr_min == xr_max:
        xr_range = [xr_min]
    else:
        xr_range = arange(xr_min,xr_max,(xr_max-xr_min)/15)
    Ub = []
    for x in xr_range:
        cap = xr_discr(xr, x, A, b, Aeq, beq)
        mcap = chq.Mobius(cap)
        if nonzero(mcap<0)[0].any():
            upbound = [[chq.max_choquet(p),p] for p in chq.cap_dnf_t(mcap,nonzero(mcap<0)[0][0],cmatr = zeros((chq.dim,chq.dim),dtype=int),result = [])]
            [q.extend([chq.Choquet(array(q[0]),q[1]) -x ,x]) for q in upbound]
            Ub.extend(upbound)
        else:
            upbound = [chq.max_choquet(cap),cap]
            upbound.extend([chq.Choquet(array(upbound[0]),upbound[1]) -x ,x])
            Ub.append(upbound)     
    return Ub,xr_min,xr_max

def caps_mm(xub_arr,xr,A,b,Aeq,beq):
    print "IN CAPS_MM"
    t0 = time()
    print len(xub_arr)
    c = array([ravel(cvx.solvers.lp(-chq.MobiusB(p)+chq.MobiusB(xr), A, b, Aeq, beq,'glpk')['x']) for p in xub_arr]).round(9)
    print "LPS", time()-t0
    c = unique(c.view([('',c.dtype)]*c.shape[1])).view(c.dtype).reshape(-1,c.shape[1])
    print "CAPS_MM", time()-t0
    return c
            
def main():
    t0 = time()
    A,b = gen_inequalities(chq.dim)
    Aeq,beq = gen_equalities(chq.dim,Shapley)
#    xr = array([1./3,1./3,1./3])
#    xr = array([0.25,0.25,0.25,0.25])
#    xr = array([0.2,0.2,0.2,0.2,0.2])
    xr = array([1./6,1./6,1./6,1./6,1./6,1./6])
    R = 0
    Umm = []
    for i in range(maxiter):
        print "XRSCAN"
        Ub,xr_min,xr_max = xr_scan(xr, A, b, Aeq, beq)
        # fix for speed
        xub_arr =  [p[0] for p in Ub if p[2] >= R]
        if not xub_arr:
            print "STOP"
            print "solution", xr
            print "distance", R
            print "active caps", Umm
            break
    #    print "XUB",xub_arr
    #    print "CAPS_MM"
    #    print caps_mm(xub_arr,xr,A,b,Aeq,beq)
        "UMM EXTEND"
        umm_bef = len(Umm)
        Umm.extend([[chq.max_choquet(p),p] for p in caps_mm(xub_arr,xr,A,b,Aeq,beq) if tuple(p) not in [tuple(q[1]) for q in Umm]])
        if len(Umm) == umm_bef:
            print "STOP NO_ADD"
            print "solution", xr
            print "distance", R
            print "active caps", Umm
            break
        print "Umm"
        print array([p[1] for p in Umm])
        print "MMAX"
        xr,R = chq.solve_mmax(Umm)
    print "solution", xr
    print "distance", R
#    print "active caps", Umm
#    print max([p[2] for p in Ub])
#    print [p[1] for p in Ub if p[2]>R]
    gap = max([p[2] for p in Ub]) - R
    print "Final and GAP", i, xr,R,gap,gap/R*100
    print t0-time() 
#    x = arange(0,1.0,0.01)
#    z = [[(-cvx.solvers.lp(chq.MobiusB(xr)-chq.MobiusB([p,q,1-p-q]),A,b,Aeq,beq,'glpk')['primal objective'],[p,q,1-p-q]) for q in x if p+q <= 1] for p in x]
#    print "LALALALA", max([p for p in chain(*z)])
#    [p for p in chain(*z) if p[0] >= R-0.005]
#    print z
#    print min(z)
#    plt.hlines(R,xr_min,xr_max)
#    plt.plot([p[3] for p in Ub],[p[2] for p in Ub],'bo')
#    plt.show()

def checker():
    A,b = gen_inequalities(chq.dim)
    Aeq,beq = gen_equalities(chq.dim,Shapley)
    xr = [0.17094965160498171, 0.18331597234367891, 0.10839328969331642, 0.17911369545266603, 0.17911369545268821, 0.17911369545266886]
    R = 0.180466679643
    for i in range(1000):
        xpick = [rnd.random() for i in range(chq.dim)]
        xpick = xpick/sum(xpick)
        z = -cvx.solvers.lp(chq.MobiusB(xr)-chq.MobiusB(xpick),A,b,Aeq,beq,'glpk')['primal objective']
        if z >= 0.16:
            print z


#import matplotlib.pyplot as plt
main()
#checker()





