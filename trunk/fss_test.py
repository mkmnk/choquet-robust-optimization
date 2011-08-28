import sys
sys.path.append('/home/gb/Downloads/science/algencan-2.3.7/bin/py')
import pdb
from numpy import *
from scipy import *
from scipy import sparse
from itertools import *
from openopt import LP,NLP
from math import *
import random as rnd
from polyhedron import Vrep, Hrep
import cvxopt as cvx
import operator
import Choquet_toolpack as chq
from capacity_parameters import *
from time import time
#from enthought.mayavi.mlab import *

set_printoptions(edgeitems='Inf',linewidth=200,suppress=True)
cvx.printing.options['width'] = -1
cvx.solvers.options['LPX_K_MSGLEV'] = 0

n = 2**chq.dim
v = zeros(n)

Shapley = {'Sh_values': [(0b001,0.2),(0b010,0.3),(0b100,0.1)], 
           'Sh_order': [],#[(0b001,0b010),(0b010,0b100)],#(0b1000,0b0100)], 
           'Sh_equal': [],#[(0b0010,0b0001)],
           'Sh_delta': 0.05}

Int_index = {'ii_values': [], 
           'ii_order': [],
           'ii_positive': [],#[(0b0001,0b0010),(0b0001,0b0100),(0b0001,0b1000),(0b0010,0b0100),(0b0010,0b1000),(0b1000,0b0100)],  
           'ii_negative': [],
           'ii_equal': [],
           'ii_delta': 0.00}
Necessity = [] 

A,b = gen_inequalities(chq.dim,Shapley,Int_index,convex=0)
Aeq,beq = gen_equalities(chq.dim,Shapley['Sh_values'],Int_index,Necessity, k_additive=2)
A = matrix(cvx.matrix([A,Aeq,-Aeq]))
b = matrix(cvx.matrix([b,beq,-beq]))
print hstack([b, -A])
print shape(A)
# raw_input()
# print "to hrep"
# a = time()
# p = Hrep(A, b)
# print "to vrep"
# Vm = array(p.generators)
# print Vm
# print a - time()
# print shape(Vm)
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


def find_cap(cap_array):
    Budg = 1  
    x0 = array([1./shape(cap_array)[0] for i in range(shape(cap_array)[0])])
#    x0 = array([1., 0, 0, 0, 0, 0, 0 ])
    print x0
    Aeq = ones(shape(cap_array)[0])
    beq = 1
    A = diag(-ones(shape(cap_array)[0]))
    b = zeros(shape(cap_array)[0])
    vmax = array([chq.Choquet(array(chq.max_choquet(v,Budg)),v) for v in cap_array])
    
    def fprime(x,cap_array,Budg):
        cap_int = dot(cap_array.T,x)
        x_int = array(chq.max_choquet(cap_int,Budg))
        ch_xr = array([chq.Choquet(x_int,v) for v in cap_array])
        ch_diffs = ch_xr - vmax
        return ch_diffs

    def objfunc(x,cap_array,Budg):
        cap_int = dot(cap_array.T,x)
        x_int = array(chq.max_choquet(cap_int,Budg))
        ch_xr = array([chq.Choquet(x_int,v) for v in cap_array])
        ch_diffs = ch_xr - vmax
#        print  dot(x, ch_diffs)
        return dot(x, ch_diffs)
    
    objfunc_w = lambda x: objfunc(x,cap_array,Budg)
    fprime_w = lambda x: fprime(x,cap_array,Budg)
    
    p = NLP(objfunc_w, x0, df=fprime_w,  A=A,  b=b,  Aeq=Aeq,  beq=beq, iprint = 25, maxIter = 10000, maxFunEvals = 1e7, diffint=1e-15, xtol=1e-11,ftol=1e-11,gtol=1e-11)
    r = p.solve('algencan')
    print r.xf
    cap_mix = dot(cap_array.T, array(r.xf))
    return cap_mix

# cap_mix = find_cap(Vm)
# x_cap_mix = chq.max_choquet(cap_mix)
# print "DIFFS ROBUST CAP"
# print array([chq.Choquet(array(chq.max_choquet(v,1)),v) - chq.Choquet(array(x_cap_mix),v) for v in Vm])
# print "ROBUST CAP"
# print cap_mix

#print "Time", time() - t0

# feasibilty test
# xr = [1./chq.dim for i in range(chq.dim)]
# xr_max = -cvx.solvers.lp(-chq.MobiusB(xr),A,b,Aeq,beq,'glpk')['primal objective']
# xr_min = cvx.solvers.lp(chq.MobiusB(xr),A,b,Aeq,beq,'glpk')['primal objective']
# print xr_max,xr_min




# for v in Vm:                 
#     d = dot(Ac,v)            
#     print d                  
#     if shape(nonzero(d<-0.00000001)[0])[1]:
#         print "NONCONVEX"    
#         # mcap = chq.Mobius(v)
#         # print chq.cap_dnf_t(mcap,nonzero(mcap<0)[0][0],cmatr = zeros((chq.dim,chq.dim),dtype=int),result =[])

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
    Z = zeros([n,n])
    for i in range(n):
        for j in range(i+1):
            if (j & i == j):
                Z[j,i] = 1
    return Z

def ncmmax(Vm):
    sols = []
    n = int(log2(shape(Vm)[1]))
    Ac = convexity(int(pow(2,n)))
    Ac = matrix(cvx.matrix(Ac))
    base = [int(pow(2,i)) for i in range(n)]
    Vm_conv = []
    Vm_nonconv = []
    for i in range(shape(Vm)[0]):
        if shape(nonzero(dot(Ac,Vm[i,:])<-0.00000001)[0])[1]:
            Vm_nonconv.append(Vm[i,:])
        else:
            Vm_conv.append(Vm[i,:])
    Uc_conv = [[chq.max_choquet(p),p] for p in Vm_conv]
    Z = zeta_matrix(int(pow(2,n)))
    for perm in permutations(base):
        print perm
        Vc = zeros(shape(Vm_nonconv))
        for j in range(n):
            if(j == 0):
                elprev = 0
            else:
                elprev = reduce(operator.or_,perm[0:j])
                el = elprev | perm[j]
                Vc[:,el] = Vm[:,el] - Vm[:,elprev]
        Vc = dot(Vc,Z)
        Uc_nonconv = [[chq.max_choquet(p),p] for p in Vc]
        Uc_nonconv.extend(Uc_conv)
        sols.append(chq.solve_mmax(Uc_nonconv))        
    return sols
#sols = ncmmax(Vm)

def convert2kadd(A,b,k=2):
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

A, b, bas = convert2kadd(A,b)
bas.insert(0,0)
A_red = hstack([b, -A])
print A_red
print shape(A_red)

# p = Hrep(A, b)
# print "to vrep"
# Vm = array(p.generators)
# print Vm
# print shape(Vm)



###### test results
##### n=4
### 1-0.2 (60, 16)
# [0.25000000000195693,
#    0.24999999999933609,
#    0.24999999999933989,
#    0.24999999999936712],
#   0.25374453730437979
### 1-0.2, 2-0.3 (68, 16)
# ([0.23200052734690277,
#    0.36715555927081928,
#    0.20042195669113913,
#    0.20042195669113888],
#   0.18023785440450746)
### 1-0.2, 2-0.3, 3-0.1 (72, 16)
# ([0.26759760202583915,
#    0.19720710515208287,
#    0.26759769149070389,
#    0.26759760133137411],
#   0.080196385759149863)

##### n=5
### 1-0.2 (140)
# [0.19999999999993406,
#    0.20000000000001653,
#    0.2000000000000163,
#    0.20000000000001641,
#    0.20000000000001664],
#   0.33018962065833857)]
###
### 1-0.2, 2-0.3 (228)
# [0.15389294589885175,
#    0.15389294589746672,
#    0.15389294590364666,
#    0.3844282163963883,
#    0.15389294590364649],
#   0.26581136651920195)]
