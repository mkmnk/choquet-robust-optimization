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
#import psyco
#psyco.full()
set_printoptions(edgeitems='Inf',linewidth=90)
#cvx.printing.options['dformat'] = '%.3f'
cvx.printing.options['width'] = -1
cvx.solvers.options['LPX_K_MSGLEV'] = 0
#### Parameters
maxiter = 10
GRID=10
#Shapley = [(0b00001,0.2),(0b00010,0.3),(0b00100,0.1),(0b01000,0.2)]
#Shapley = [(0b00001,0.3),(0b00010,0.3),(0b00100,0.1)]
#Shapley = [(0b00001,0.5),(0b00010,0.3),(0b00100,0.1)]
#Shapley = [(0b00001,0.2),(0b00010,0.3),(0b00100,0.1)]
#Shapley = [(0b00001,0.7),(0b00010,0.1)]
#Shapley = [(0b00001,0.2),(0b00010,0.3)]
#Shapley = [(0b00100,0.8)]
#Shapley = [(0b00010,0.4)]
#Shapley = [(0b00001,0.2)]
#Shapley = []

# example thesis
# NODE 4
#00001 - fw
#00010 - part
#00100 - NIDS
#01000 - IPS
#10000 - Encr
#Shapley = {'Sh_values': [(0b00001,0.4)], 
#           'Sh_order': [(0b00001,0b00010),(0b00010,0b10000),(0b10000,0b00100)], 
#           'Sh_equal': [(0b00100,0b01000)],
#           'Sh_delta': 0.05}
#
#Int_index = {'ii_values': [], 
#           'ii_order': [((0b00001,0b01000),(0b00001,0b00100))],
#           'ii_positive': [(0b00001,0b00100)], 
#           'ii_equal': [],
#           'ii_delta': 0.05}
#Necessity = []
## NODE 10
#Shapley = {'Sh_values': [], 
#           'Sh_order': [], 
#           'Sh_equal': [],
#           'Sh_delta': 0.05}
#
#Int_index = {'ii_values': [], 
#           'ii_order': [],
#           'ii_positive': [(0b0010,0b1000)], 
#           'ii_equal': [],
#           'ii_delta': 0.05}
#Necessity = [0b0001] 
#
#
## NODE 17
#Shapley = {'Sh_values': [(0b100,0.45)], 
#           'Sh_order': [], 
#           'Sh_equal': [],
#           'Sh_delta': 0.05}
#
#Int_index = {'ii_values': [], 
#           'ii_order': [],
#           'ii_positive': [], 
#           'ii_equal': [],
#           'ii_delta': 0.05}
#Necessity = [0b100, 0b010] 
# NODE 21, 22
#Shapley = {'Sh_values': [], 
#           'Sh_order': [(0b10,0b01)], 
#           'Sh_equal': [],
#           'Sh_delta': 0.05}
#
#Int_index = {'ii_values': [], 
#           'ii_order': [],
#           'ii_positive': [], 
#           'ii_equal': [],
#           'ii_delta': 0.05}
#Necessity = [0b10] 

# NODE 23
#00001 - sec policy
#00010 - acc cont
#00100 - user mgmt
#01000 - pol compl
#10000 - str auth
#Shapley = {'Sh_values': [], 
#           'Sh_order': [(0b00010,0b00100),(0b00001,0b01000),(0b01000,0b10000)], 
#           'Sh_equal': [(0b00100,0b00001)],
#           'Sh_delta': 0.05}
#
#Int_index = {'ii_values': [], 
#           'ii_order': [],
#           'ii_positive': [(0b00100,0b01000)], 
#           'ii_equal': [],
#           'ii_delta': 0.05}
#Necessity = [0b00010]

## NODE 30
#Shapley = {'Sh_values': [(0b010,0.5)], 
#           'Sh_order': [(0b010,0b100),(0b100,0b001)], 
#           'Sh_equal': [],
#           'Sh_delta': 0.05}
#
#Int_index = {'ii_values': [], 
#           'ii_order': [],
#           'ii_positive': [], 
#           'ii_equal': [],
#           'ii_delta': 0.05}
#Shapley = {'Sh_values': [(0b00001,0.4)], 
#           'Sh_order': [(0b00001,0b00010),(0b00010,0b10000),(0b10000,0b00100)], 
#           'Sh_equal': [(0b00100,0b01000)],
#           'Sh_delta': 0.05}
#
#Int_index = {'ii_values': [], 
#           'ii_order': [((0b00001,0b01000),(0b00001,0b00100))],
#           'ii_positive': [(0b00001,0b00100)], 
#           'ii_equal': [],
#           'ii_delta': 0.05}
#Necessity = []
#Necessity = [] 

### NODE 3 Approximated F(B) for 4 and 10
### 4net 10ser 15phy 16ws
#Shapley = {'Sh_values': [(0b0010,0.4)], 
#           'Sh_order': [(0b0001,0b1000),(0b1000,0b0100)], 
#           'Sh_equal': [(0b0010,0b0001)],
#           'Sh_delta': 0.05}
#
#Int_index = {'ii_values': [], 
#           'ii_order': [],
#           'ii_positive': [], 
#           'ii_equal': [],
#           'ii_delta': 0.05}
#Necessity = [0b0001, 0b0010] 

## NODE 22 Approximated F(B) for 23 and NODE 21 Approximated F(B) for 22 and 30
## 23Acc 29Aud  ## 22TecSt 30HR
#Shapley = {'Sh_values': [], 
#           'Sh_order': [(0b01,0b10)], 
#           'Sh_equal': [],
#           'Sh_delta': 0.05}
#
#Int_index = {'ii_values': [], 
#           'ii_order': [],
#           'ii_positive': [], 
#           'ii_equal': [],
#           'ii_delta': 0.05}
#Necessity = [0b01]

# rcap test
Shapley = {'Sh_values': [(0b001,0.2)], 
           'Sh_order': [], 
           'Sh_equal': [],
           'Sh_delta': 0.05}

Int_index = {'ii_values': [], 
           'ii_order': [],
           'ii_positive': [], 
           'ii_equal': [],
           'ii_delta': 0.05}
Necessity = [] 


#####
#
######MAIN



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
        xr_range = arange(xr_min,xr_max,(xr_max-xr_min)/GRID)
    Ub = []
    for x in xr_range:
        cap = xr_discr(xr, x, A, b, Aeq, beq)
        mcap = chq.Mobius(cap)
#        mcap = array([ 0. ,  0.6,  0.6,  0.6, -0.3, -0.3, -0.3,  0.1])
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
    t0 = time()
    A,b = gen_inequalities(chq.dim,Shapley,Int_index)
    Aeq,beq = gen_equalities(chq.dim,Shapley['Sh_values'],Int_index,Necessity)
#    print A, Aeq
#    raw_input("123")
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
    #    print "XUB",xub_arr
    #    print "CAPS_MM"
    #    print caps_mm(xub_arr,xr,A,b,Aeq,beq)
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
#    print "active caps", Umm
#    print max([p[2] for p in Ub])
#    print [p[1] for p in Ub if p[2]>R]
    gap = max([p[2] for p in Ub]) - R
    print "Final and GAP", i, xr,R,gap,gap/R*100
    print "Time", time()-t0 
    return array(xr)
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

xr = main()


#import matplotlib.pyplot as plt
#checker()





