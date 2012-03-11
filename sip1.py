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
from capacity_parameters import *
cvx.solvers.options['show_progress'] = False
import matplotlib.pyplot as plt
#psyco.full()
set_printoptions(edgeitems='Inf',linewidth=90)
#cvx.printing.options['dformat'] = '%.3f'
cvx.printing.options['width'] = -1
cvx.solvers.options['LPX_K_MSGLEV'] = 0
#### Parameters
maxiter = 10
GRID=21

#####
#parameters
####
from node_config import *
###


def xr_discr(xr, threshold, A, b, Aeq, beq):
    """
    Input: xr - current point, threshold - value in [Cmin,Cmax] at xr, U = A*v<b,Aeq*v=beq
    calculates maximal (pointwise) capacity v^max in U, such that (threshold=C(v,xr)). Used for discretization of [Cmin, Cmax] interval
    """
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
    """
    Input: xr - current point, U = A*v<b,Aeq*v=beq, GRID(global)
    Calculates Cmin,Cmax and GRID points in between. For each point
    Calculates v^max and decomposes it, calculates all maxima
    Output: regret values and points
    """
    print "XR",xr
    xr_max = -cvx.solvers.lp(-chq.MobiusB(xr),A,b,Aeq,beq,'glpk')['primal objective']
    xr_min = cvx.solvers.lp(chq.MobiusB(xr),A,b,Aeq,beq,'glpk')['primal objective']
    if xr_min == xr_max:
        xr_range = [xr_min]
    else:
        xr_range =linspace(xr_min,xr_max,GRID)
    Ub = []
    for x in xr_range:
        cap = xr_discr(xr, x, A, b, Aeq, beq)
        mcap = chq.Mobius(cap)
        if nonzero(mcap<0)[0].any():
            """ Overcutting
            c = array([chq.max_choquet(p) for p in chq.cap_dnf_t(mcap,nonzero(mcap<0)[0][0],cmatr = zeros((chq.dim,chq.dim),dtype=int),result = [])]).round(2)
            print "before\n",c
            c = unique(c.view([('',c.dtype)]*c.shape[1])).view(c.dtype).reshape(-1,c.shape[1])
            print "after\n",c
            print mcap
            raw_input("Press Enter to continue...")
            """
            upbound = [[chq.max_choquet(p),p] for p in chq.cap_dnf_t(mcap,nonzero(mcap<0)[0][0],cmatr = zeros((chq.dim,chq.dim),dtype=int),result = [])]
            [q.extend([chq.Choquet(array(q[0]),q[1]) -x ,x]) for q in upbound]
            Ub.extend(upbound)
        else:
            upbound = [chq.max_choquet(cap),cap]
            upbound.extend([chq.Choquet(array(upbound[0]),upbound[1]) -x ,x])
            Ub.append(upbound)     
    return Ub,xr_min,xr_max



def caps_mm(xub_arr,xr,A,b,Aeq,beq,R):
    """
    Input: upper bound array, current point xr, U=A*v<b,Aeq*v=beq, current regret R
    Calc max regret over U at point from xubb_arr
    Output: unique capacities which produce regret > R at points p from xub_arr
    """
    print "IN CAPS_MM"
    t0 = time()
    print len(xub_arr)
    Sols = [cvx.solvers.lp(-chq.MobiusB(p)+chq.MobiusB(xr), A, b, Aeq, beq,'glpk') for p in xub_arr]   
    c = array([ravel(p['x']) for p in Sols if p['primal objective'] < -R]).round(13)
    if len(c) == 0:
        print "NO CAPSMM"
        return c
    print "LPS", time()-t0
    c = unique(c.view([('',c.dtype)]*c.shape[1])).view(c.dtype).reshape(-1,c.shape[1])
    print "CAPS_MM", time()-t0
    return c
            


def main():
    """
    Input:none
    Generate U=A*v<b,Aeq*v<beq
    do a Cmin to Cmax scan at xr, find upper bound candidate points
    find capacities which exceed upper bound in these points, add to Umm
    solve mmax for Umm
    Output: solution, regrets for all v in Umm, plot
    """
    t0 = time()
    A,b = gen_inequalities(chq.dim,Shapley,Int_index,convex=1)
    Aeq,beq = gen_equalities(chq.dim,Shapley['Sh_values'],Int_index,Necessity)
    xr = [1./chq.dim for i in range(chq.dim)]
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
        "UMM EXTEND"
        umm_bef = len(Umm)
        Umm.extend([[chq.max_choquet(p),p] for p in caps_mm(xub_arr,xr,A,b,Aeq,beq,R) if tuple(p) not in [tuple(q[1]) for q in Umm]])
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
    print "active distance"
    for i in [(p[1], chq.Choquet(array(p[0]),p[1]) - chq.Choquet(array(xr),p[1])) for p in Umm]:
        print i
    gap = max([p[2] for p in Ub]) - R
    gap_arr =  [p for p in Ub if p[2] >= R]
    print "Final and GAP", i, xr,R,gap,gap/R*100
    print "Time", time()-t0
    plt.hlines(R,xr_min,xr_max)
    plt.plot([p[3] for p in Ub],[p[2] for p in Ub],'bo')
    plt.show() 
    return array(xr),Ub


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

xr,Ub = main()


 

