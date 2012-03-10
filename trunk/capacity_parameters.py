'''
Created on 11 Mar 2011

@author: gb
'''

from numpy import *
from itertools import *
from math import factorial
import random as rnd
import cvxopt as cvx
import Choquet_toolpack as chq

def convexity(m):
    """
    Input: m = 2^dim. Output: matrix with v(AuB) + v(AnB) - v(A) - v(B) for all A,B \in N
    """
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

def monotonicity(m):
    """
    Input: m = 2^dim. Output: matrix with v(A) - v(B) for all B \in A \in N
    """
    x = []
    I = []
    J = []
    j = 0
    for supset in range(m):
        for subset in range(supset):
            if (supset & subset == subset) & (chq.bitCount(supset) - chq.bitCount(subset) == 1):
                x.extend([1,-1])
                I.extend([j,j])
                J.extend([supset,subset])
                j = j+1
    A = cvx.spmatrix(x,I,J)
    return A

def k_additivity(m,k):
    """
    Input: m=2^dim, k \in [1,dim] Output: mobius expressions for all A \in N, |A| > k
    """
    x = []
    I = []
    J = []
    j = 0
    for supset in range(m):
        if (chq.bitCount(supset) > k):
            for subset in range(supset):
                if (supset & subset == subset):
                    x.extend([pow(-1,chq.bitCount(supset-subset))])
                    I.extend([j])
                    J.extend([subset])
            x.extend([1])
            I.extend([j])
            J.extend([supset])
            j = j+1
    A = cvx.spmatrix(x,I,J)
    return A


def limits(m):
    """
   Input: m=2^dim, Output: matrix of ones on a diagonal 
    """
    x = ones(m-2)
    I = range(m-2)
    J = range(1,m-1)   
    A = cvx.spmatrix(x,I,J,(m-2,m))
    return A

def equalities(m):
    """
    Input:m=2^dim Output: v(0)=0, v(1)=1 of a correct size
    """
    A = cvx.spmatrix([1,1],[0,1],[0,m-1],(2,m))
    return A

def shapley(m,el):
    """
    Input: m=2^dim, el \in [0,dim]. Output: Shapley coefficient vector S for a particular criterion el. <S,v> = sh(el)
    """
    x = []
    I = []
    J = []
    for i in range(m):        
        if i & el != el:
            coeff = float(factorial(chq.bitCount(m-1) - chq.bitCount(i) - 1))*factorial(chq.bitCount(i))/factorial(chq.bitCount(m-1))
            x.extend([coeff,-coeff])
            I.extend([0,0])   # for some reason it wouldn't take I=zeros(m)
            J.extend([i | el, i])
    A = cvx.spmatrix(x,I,J)
    return A  

def gen_inequalities(dim, Shapl_dict = {}, II_dict = {}, convex=0):
    """
    Input: dim, Shapley values structure, II values structure, convexity flag
    calls the appropriate procedures
    Output: Matrix and a column, Av<b
    """
    if convex:
        A = cvx.sparse([-convexity(2**dim),-limits(2**dim)])
        b = cvx.matrix(0,(A.size[0],1),"d")
    else:
        A = cvx.sparse([-monotonicity(2**dim),-limits(2**dim)])
        b = cvx.matrix(0,(A.size[0],1),"d")
    if Shapl_dict['Sh_order']:
        Sh_diffs = - cvx.matrix([shapley(2**dim,p[0]) - shapley(2**dim,p[1]) for p in Shapl_dict['Sh_order']])
        b_sh_diffs = - cvx.matrix(Shapl_dict['Sh_delta'],(Sh_diffs.size[0],1),"d")
        A = cvx.sparse([A, Sh_diffs])
        b = cvx.matrix([b,b_sh_diffs])
    if Shapl_dict['Sh_equal']:
        Sh_eq1 = cvx.matrix([shapley(2**dim,p[0]) - shapley(2**dim,p[1]) for p in Shapl_dict['Sh_equal']])
        b_sh_eq1 = cvx.matrix(Shapl_dict['Sh_delta'],(Sh_eq1.size[0],1),"d")
        # -Sh_eq2 < delta
        Sh_eq2 = - cvx.matrix([shapley(2**dim,p[0]) - shapley(2**dim,p[1]) for p in Shapl_dict['Sh_equal']])
        b_sh_eq2 = cvx.matrix(Shapl_dict['Sh_delta'],(Sh_eq2.size[0],1),"d")
        A = cvx.sparse([A, Sh_eq1, Sh_eq2])
        b = cvx.matrix([b, b_sh_eq1, b_sh_eq2])
    if II_dict['ii_values']:
        II_vals = - cvx.matrix([int_index(2**dim,p[0]) for p in II_dict['ii_values']])
        b_ii_vals = - cvx.matrix([p[1] for p in II_dict['ii_values']])
        A = cvx.sparse([A, II_vals])
        b = cvx.matrix([b,b_ii_vals])
    if II_dict['ii_order']:
        II_diffs = - cvx.matrix([int_index(2**dim,p[0]) - int_index(2**dim,p[1]) for p in II_dict['ii_order']])
        b_ii_diffs = - cvx.matrix(II_dict['ii_delta'],(II_diffs.size[0],1),"d")
        A = cvx.sparse([A, II_diffs])
        b = cvx.matrix([b,b_ii_diffs])
    if II_dict['ii_positive']:
        II_positive = -cvx.matrix([int_index(2**dim,p) for p in II_dict['ii_positive']])
        b_ii_positive = - cvx.matrix(II_dict['ii_delta'],(II_positive.size[0],1),"d")
        A = cvx.sparse([A, II_positive])
        b = cvx.matrix([b, b_ii_positive])
    if II_dict['ii_negative']:
        II_negative = cvx.matrix([int_index(2**dim,p) for p in II_dict['ii_negative']])
        b_ii_negative = - cvx.matrix(II_dict['ii_delta'],(II_negative.size[0],1),"d")
        A = cvx.sparse([A, II_negative])
        b = cvx.matrix([b, b_ii_negative])
    return A,b

def gen_equalities(dim,Shpl,II = [], Necessity = [], Sufficiency = [], k_additive=0 ):
    """
    Input: dim, Shapley values structure, II values structure, Necessity criteria, Sufficiency criteria,k-additivity
    calls the appropriate procedures
    Output: Matrix and a column, Aeq*v=beq
    """
    Alist = [equalities(2**dim)]
    Alist.extend([shapley(2**dim,p[0]) for p in Shpl])
    Aeq = cvx.sparse(Alist)
    blist = [0.,1]
    blist.extend([p[1] for p in Shpl])
    beq = cvx.matrix(blist)
    if Necessity:
        Aeq = cvx.sparse([Aeq, necessity(2**dim,Necessity)])
        beq = cvx.matrix([beq, cvx.matrix(0.,(size(Aeq)[0]-size(beq)[0],1))])
    if Sufficiency:
        Aeq = cvx.sparse([Aeq, sufficiency(2**dim,Sufficiency)])
        beq = cvx.matrix([beq, cvx.matrix(1.,(size(Aeq)[0]-size(beq)[0],1))])
    if k_additive:
        Aeq = cvx.sparse([Aeq, k_additivity(2**dim,k_additive)])
        beq = cvx.matrix([beq, cvx.matrix(0.,(size(Aeq)[0]-size(beq)[0],1))])
    return Aeq,beq

def int_index(m,el_pair):
    """
    Input: m=2^dim, el_pair \in [0,dim]^2. Output: II coefficient vector II for a particular criterion el. <II,v> = II(el1,el2)
    """
    x = []
    I = []
    J = []
    for i in range(m):        
        if (i & el_pair[0] != el_pair[0]) and (i & el_pair[1] != el_pair[1]):
            coeff = float(factorial(chq.bitCount(m-1) - chq.bitCount(i) - 2))*factorial(chq.bitCount(i))/factorial(chq.bitCount(m-1)-1)
            x.extend([coeff,-coeff,-coeff,coeff])
            I.extend([0,0,0,0])   # for some reason it wouldn't take I=zeros(m)
            J.extend([i | el_pair[0] | el_pair[1], i | el_pair[0], i | el_pair[1], i])
    A = cvx.spmatrix(x,I,J)
    return A

def necessity(m,els):
    """
    Input: m=2^dim, els - list of neccessary subsets. Output: matrix with ones for A not including any of els. 1 set per row
    """
    x = []
    I = 0
    J = []
    for i in range(1,m):        
        if ([i & p for p in els] != els):
            x.extend([1])
            I = I+1  
            J.extend([i]) 
    A = cvx.spmatrix(x,range(I),J,(I,m))
    return A

def sufficiency(m,els):
    """
    Input: m=2^dim, els - list of sufficient subsets. Output: matrix with ones for A including any of els. 1 set per row
    """
    x = []
    I = 0
    J = []
    for i in range(1,m):        
        for p in els:
            if (i & p == p):
                x.extend([1])
                I = I+1  
                J.extend([i]) 
    A = cvx.spmatrix(x,range(I),J,(I,m))
    return A
