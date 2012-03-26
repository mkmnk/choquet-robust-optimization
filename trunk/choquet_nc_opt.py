import sys
from numpy import *
from scipy import *
from scipy import sparse
from itertools import *
from openopt import LP,NSP
from math import *
import random as rnd
from polyhedron import Vrep, Hrep
import cvxopt as cvx
import operator
import Choquet_toolpack as chq
from capacity_parameters import *
from time import time
import Choquet_toolpack
sys.path.append('/home/gb/Downloads/science/algencan-2.3.7/bin/py')
#from enthought.mayavi.mlab import *

set_printoptions(edgeitems='Inf',linewidth=200,suppress=True)
# cvx.printing.options['width'] = -1
# cvx.solvers.options['LPX_K_MSGLEV'] = 0

# n = 2**chq.dim
# v = zeros(n)

#from node_config import *
# Vx = [(chq.max_choquet(p),p) for p in Vm]
# print Vx
# print "MMAX"
# xr = chq.solve_mmax(Vx)
# print xr
# #print Vx
# print array([chq.Choquet(array(xr[0]),v) - chq.Choquet(array(chq.max_choquet(v,1)),v) for v in Vm])
# cap_r = array([0.000000,0.383333,0.333333,0.716667,0.283333,0.666667,0.616667,1])
# x_cap_r = chq.max_choquet(cap_r)
# print "Entropy cap diff"
# print array([chq.Choquet(array(chq.max_choquet(v,1)),v) - chq.Choquet(array(x_cap_r),v) for v in Vm])


# cap_mix = solve_dual(Vm)
# x_cap_mix = chq.max_choquet(cap_mix)
# print "DIFFS ROBUST CAP"
# print array([chq.Choquet(array(chq.max_choquet(v,1)),v) - chq.Choquet(array(x_cap_mix),v) for v in Vm])
# print "ROBUST CAP"
# print cap_mix


def zeta(cap):
    """
    Inverse to mobius
    """
    nu = zeros(len(cap))
    i = 0
    j = 0
    for i in range(1,len(cap)):
        for j in range(0,i+1):
            if i & j == j:
                nu[i] += cap[j]
    return nu

def zeta_matrix(n):
    """
    Outputs a matrix which multiplied by a matrix of Mobius coeffs gives a matrix of capacities
    """
    Z = zeros([n,n])
    for i in range(n):
        for j in range(i+1):
            if (j & i == j):
                Z[j,i] = 1
    return Z

def nc_primal(Vm):
    """
    Nonconvex minimax - for a matrix of capacities finds all nonconvex rows, and for all permutations generates necessity measures for them. Then solves mmax for the resulting matrix
    """
    sols = []
    n = int(log2(shape(Vm)[1]))
    Ac = convexity(int(pow(2,n)))
    Ac = matrix(cvx.matrix(Ac))
    base = [int(pow(2,i)) for i in range(n)]
    Vm_conv = []
    Vm_nonconv = []
    for i in range(shape(Vm)[0]):                                          # test for 2-monotonicity, separate into 2 arrays
        if shape(nonzero(dot(Ac,Vm[i,:])<-0.00000001)[0])[1]:
            Vm_nonconv.append(Vm[i,:])
        else:
            Vm_conv.append(Vm[i,:])
    Vm_conv=array(Vm_conv)
    Vm_nonconv=array(Vm_nonconv)
    Uc_conv = [[chq.Choquet(array(chq.max_choquet(p)),p),p] for p in Vm_conv]
    Uc_nonconv = [[chq.Choquet(array(chq.max_choquet_glob(p)),p),p] for p in Vm_nonconv]
    Z = zeta_matrix(int(pow(2,n)))
    for perm in permutations(base):
        print perm
        Uc_necmeas = [[p[0],nec_measure(p[1],perm,Z)] for p in Uc_nonconv]
        Uc_necmeas.extend(Uc_conv)
        Gmax = array([i[0] for i in Uc_necmeas])
        Vconv = array([i[1] for i in Uc_necmeas])
        rcap,rmix = Choquet_toolpack.solve_dual(Vconv,Gmax)
        sol = list(chq.solve_mmax_wval(Uc_necmeas))
        sol.extend([chq.max_choquet(rcap),rcap, perm])
        sols.append(sol)
    return sols
    

def convert2kadd(A,b,k=2):
    """
    For a k-additive capacity v, all v(A) for |A|>k can be expressed through v(B), B in A, |B|<=k
    In other words certain coordinates of v are linearly dependent
    We transform A so that A*v<b becomes Ai * v(<=k) < b
    This reduces the dimension of A (and accordingly v) significantly and eases extr. point search
    for 2-add capacities (this implementation) v(A) = sum(v(i,j)) - (|A|-2)sum(v(i)) i,j \in A
    All columns to be multiplied by v(A), |A|<=2 are updated according to the equality above
    bas is array of linearly independent column indexes
    """
    n = shape(A)[1]
    bas = []
    for i in range(1,n):
        card_i = chq.bitCount(i)
        if (card_i <= k):
            for j in range(i,n):
                card_j = chq.bitCount(j)
                if (i & j == i) & (card_j > k):
                    if (card_i == 2):
                        A[:,i] = A[:,i] + A[:,j]
                    elif (card_i == 1):
                        A[:,i] = A[:,i] - (card_j - 2)*A[:,j]
                    else:
                        print "something wrong here", i, card_i
                        break
            bas.append(i)
        else:
            A[:,i] = 0
    print A
    print bas
    A = A[:,bas]
    b = b[unique(ravel(nonzero(A)[0]))]
    A = A[unique(ravel(nonzero(A)[0])),:]
    return A,b,bas
                

def convertnadd(Vm, bas,k=2):
    """
    Extend A to include factors (i.e. columns) for linearly dependent coordinates of v.
    Inverse to convert2kadd
    This implelentation only supports k=2
    """
    n = 2**chq.dim # fix later to derive as solution of n(n+1) = 2*len(Vm)
    Vm_new = zeros([shape(Vm)[0],n])
    Vm_new[:,bas] = Vm
    for i in range(n+1):
        card_i = chq.bitCount(i)
        if (card_i > k):
            for j in range(i):
                card_j = chq.bitCount(j)
                if (j & i == j):
                    if (card_j == 2):
                        Vm_new[:,i] = Vm_new[:,i] + Vm_new[:,j]
                    elif (card_j == 1):
                        Vm_new[:,i] = Vm_new[:,i] - (card_i - 2)*Vm_new[:,j]
    return Vm_new

def nec_measure(capacity,perm,Z):
    """
    Generate necessity measure (or beta-measure) for a given permutation
    calculates Mobius and then does zeta
    """
    necmeas = zeros(len(capacity))
    for j in range(len(perm)):
        if(j == 0):
            elprev = 0
        else:
            elprev = reduce(operator.or_,perm[0:j])
        el = elprev | perm[j]
        necmeas[el] = capacity[el] - capacity[elprev]
    return dot(necmeas,Z)


def nc_dual(Vm,dm_info,Budg):
    """
    dual solving for nonconvex capacities
    actually a wrapper around solve_dual (calls it for each permutation)
    """
    sols = []
    dim = len(dm_info['criteria_functions'])
    Ac = convexity(int(pow(2,dim)))
    Ac = matrix(cvx.matrix(Ac))
    base = [int(pow(2,i)) for i in range(dim)]
    Vm_conv = []
    Vm_nonconv = []
    for i in range(shape(Vm)[0]):                                          # test for 2-monotonicity, separate into 2 arrays
        if shape(nonzero(dot(Ac,Vm[i,:])<-0.00000001)[0])[1]:
            Vm_nonconv.append(Vm[i,:])
        else:
            Vm_conv.append(Vm[i,:])
    Vm_conv=array(Vm_conv)
    Vm_nonconv=array(Vm_nonconv)
    print "NONCONVEX VERTICES"
    print "***********************************************"
    print Vm_nonconv
    print "CONVEX VERTICES"
    print "***********************************************"
    print Vm_conv
    fx = dm_info['criteria_functions'].values()
    dfx = dm_info['criteria_fgrad'].values()
    Uc_conv = [[chq.Choquet(array(chq.max_choquet(p,fx,dfx,Budg)),p,fx),p] for p in Vm_conv] 
    print "Choquet integral VALUES for convex caps"
    print "***********************************************"
    print Uc_conv
    Uc_nonconv = [[chq.Choquet(array(chq.max_choquet_glob(p,fx,dfx,Budg)),p,fx),p] for p in Vm_nonconv] # [(glob_max,capacity)]
    print "Choquet integral VALUES for nonconvex caps"
    print "***********************************************"
    print Uc_nonconv
    Z = zeta_matrix(int(pow(2,dim)))
    # SIP part (capacity filter)
    Uc_c = []
    Uc_nc = []
    zr = Choquet_toolpack.find_centre(fx,dfx)[0]   # initial parameters, might be changed
    t = 0 
    while 1:
        if t == 0:                            # kickstart: assume a large t and go down until there are capacities for which the diff exceeds it
            t = 0.4
            cl_c = []
            cl_nc = []
            while not (cl_c or cl_nc):
                cl_c = sip_filter(Uc_conv,zr,fx,t)
                cl_nc = sip_filter(Uc_nonconv,zr,fx,t)
                t = t - 0.05
        else:
            cl_c = sip_filter(Uc_conv,zr,fx,t)
            cl_nc = sip_filter(Uc_nonconv,zr,fx,t)
        print "Threshold: ",t    
        print "Filtered: convex: ", len(cl_c), " Nonconvex: ", len(cl_nc)
        if (cl_c or cl_nc):
            Uc_c.extend(cl_c)
            Uc_nc.extend(cl_nc)
            for perm in permutations(base):     # nonconvexity part (permutations)
                print perm
                Uc_necmeas = [[p[0],nec_measure(p[1],perm,Z)] for p in Uc_nc]
                Uc_necmeas.extend(Uc_c)
                Gmax = array([i[0] for i in Uc_necmeas])
                Vconv = array([i[1] for i in Uc_necmeas])
                vr,obj_d,d_mix = Choquet_toolpack.solve_dual_sip(Uc_necmeas,fx,dfx,Budg)
                sols.append((obj_d,vr))
                # sols.append((obj_p,xr))
            sols.sort(key=operator.itemgetter(0),reverse=True)
            print sols
            zr = array(Choquet_toolpack.max_choquet(sols[0][1],fx,dfx,Budg))
        RSol_D = sols[0][1]
        OBJ_D = sols[0][0]
        break 
    Uc_c = []
    Uc_nc = []
    zr = Choquet_toolpack.find_centre(fx,dfx)[0]   # initial parameters, might be changed
    t = 0 
    sols = []
    while 1:                          # Now solve primal
        if t == 0:                            # kickstart: assume a large t and go down until there are capacities for which the diff exceeds it
            t = 0.4
            cl_c = []
            cl_nc = []
            while not (cl_c or cl_nc):
                cl_c = sip_filter(Uc_conv,zr,fx,t)
                cl_nc = sip_filter(Uc_nonconv,zr,fx,t)
                t = t - 0.05
        else:
            cl_c = sip_filter(Uc_conv,zr,fx,t)
            cl_nc = sip_filter(Uc_nonconv,zr,fx,t)
        print "Threshold: ",t    
        print "Filtered: convex: ", len(cl_c), " Nonconvex: ", len(cl_nc)
        if (cl_c or cl_nc):
            Uc_c.extend(cl_c)
            Uc_nc.extend(cl_nc)
            print 123
            for perm in permutations(base):     # nonconvexity part (permutations)
                print perm
                Uc_necmeas = [[p[0],nec_measure(p[1],perm,Z)] for p in Uc_nc]
                Uc_necmeas.extend(Uc_c)
                Gmax = array([i[0] for i in Uc_necmeas])
                Vconv = array([i[1] for i in Uc_necmeas])
                xr,obj_p = Choquet_toolpack.solve_mmax_wval(Uc_necmeas,fx,dfx,Budg)
                sols.append((obj_p,xr))
            sols.sort(key=operator.itemgetter(0),reverse=True)
            print sols
        RSol_P = sols[0][1]
        OBJ_P = sols[0][0]
        break 
    print "FOUND SOLUTION"
    print "dual robust point: ", RSol_D
    print "dual objective value: ", OBJ_D
    print "primal robust point: ", RSol_P
    print "primal objective value: ", OBJ_P
    return sols[0][1]


def sip_filter(cap_list,zr,fx,t):
    """
    Goes through the vertices and picks v_i such that max_z C(v_i,f(z)) - C(v_i,f(z^r)) > t
    If t = 0 (on first iteration) returns top 3 capacities from both lists
    Input: Uc_conv,Uc_nonconv - lists of pairs (max_val, capacity), zr - point to check, fx - functions, t - threshold
    Output: two lists of capacities filtered according to the rules above
    """

    c = [x[1] for x in [(i[0] - Choquet_toolpack.Choquet(zr,i[1],fx)-t,i) for i in cap_list] if x[0]>0] # strict inequalities!

    return c

    
# A, b, bas = convert2kadd(A,b)
# bas.insert(0,0)
# A_red = hstack([b, -A])
# print A_red
# print shape(A_red)
# Vm = genfromtxt('/home/gb/Documents/papers/FSS_LINZ/cdd/shapley4_1-02_red.mat')
# Vm = convertnadd(Vm,bas)


###########################
##### single-permutation
###########################
# n = int(log2(shape(Vm)[1]))
# Ac = convexity(int(pow(2,n)))
# Ac = matrix(cvx.matrix(Ac))
# base = [int(pow(2,i)) for i in range(n)]
# Vm_conv = []
# Vm_nonconv = []
# for i in range(shape(Vm)[0]):                                          # test for 2-monotonicity, separate into 2 arrays
#     if shape(nonzero(dot(Ac,Vm[i,:])<-0.00000001)[0])[1]:
#         Vm_nonconv.append(Vm[i,:])
#     else:
#         Vm_conv.append(Vm[i,:])
# Z = zeta_matrix(int(pow(2,n)))
# perm = (1,8,2,4)
# Vm_necmeas = [nec_measure(p,perm,Z) for p in Vm_nonconv]
# Vm_necmeas.extend(Vm_conv)
# Vm_necmeas=array(Vm_necmeas)
# print Vm_necmeas
# gm_conv = [chq.Choquet(array(chq.max_choquet(p)),p) for p in Vm_conv]
# gm_nonconv = [chq.Choquet(array(chq.max_choquet_glob(p)),p) for p in Vm_nonconv]
# gm_nonconv.extend(gm_conv)
# gm_nonconv = array(gm_nonconv)
# print gm_nonconv
# vr,rmix=solve_dual(Vm_necmeas,gm_nonconv)
# print chq.max_choquet(vr)

# Vm1 = array([rmix[i]*Vm[i,:] for i in range(len(rmix))])
# vgr = array([sum(Vm1[:,i]) for i in range(16)])

# for i in range(shape(Vm)[0]):                                          # test for 2-monotonicity, separate into 2 arrays
#     if shape(nonzero(dot(Ac,Vm[i,:])<-0.00000001)[0])[1]:
#         print "Global"
#         print chq.Choquet(array(chq.max_choquet_glob(Vm[i,:])),Vm[i,:]) - chq.Choquet(glob_max,Vm[i,:])
#         print "Robust"
#         print chq.Choquet(array(chq.max_choquet_glob(Vm[i,:])),Vm[i,:]) - chq.Choquet(loc_max,Vm[i,:])
#     else:
#         print "Global"
#         print chq.Choquet(array(chq.max_choquet(Vm[i,:])),Vm[i,:]) - chq.Choquet(glob_max,Vm[i,:])
#         print "Robust"
#         print chq.Choquet(array(chq.max_choquet(Vm[i,:])),Vm[i,:]) - chq.Choquet(loc_max,Vm[i,:])


##############################################
#### end
##############################################




###### test results
##### n=4
### 1-0.2 (60, 16) FC
# [0.067862569909373741,
#    0.31071247669672775,
#    0.31071247669689966,
#    0.31071247669699886],
#   0.34392421166929649)]
### 1-0.2, 2-0.3 (68, 16) FC
# [0.092573637457428859,
#    0.17918015578706173,
#    0.36412310337775328,
#    0.36412310337775622],
#   0.28563376714613425)]
### 1-0.2, 2-0.3, 3-0.1 (72, 16) FC
 # ([0.18778767972060673,
 #   0.27915457200152449,
 #   0.023332286749396542,
 #   0.5097254615284722],
 #  0.1391229510900707)]


##### n=5
### 1-0.2 (140) FC
 # ([0.032822087620647498,
 #   0.24179447809472968,
 #   0.241794478094886,
 #   0.24179447809486782,
 #   0.24179447809486912],
 #  0.43435182954133417)
###
### 1-0.2, 2-0.3 (228) FC
# ([0.35704201479785341,
#    0.29268492617022696,
#    0.11675768634389602,
#    0.11675768634406954,
#    0.1167576863439542],
#   0.33661625336665091)### 1-0.2, 2-0.3, 3-0.1 (648) FC
# [0.084977758369359724,
#    0.41377953465179151,
#    0.080479779281832378,
#    0.21038146384850179,
#    0.21038146384851467],
#   0.26476353066627256)]
