'''
Created on 6 Mar 2011

@author: gb
'''


from sip1 import *
set_printoptions(edgeitems='Inf',linewidth=90)
cvx.printing.options['dformat'] = '%.3f'
cvx.printing.options['width'] = -1


Bs = [(1.5, [0.5142862205790838, 0.47142755884183246, 0.51428622057908369])]
#Bs = [(1,  [0.35046115511896969, 0.29907768976206106, 0.35046115511896919])]
#Bs = [(1.1, [0.38379451461205039, 0.33241097077590137, 0.38379451461204844])]
#Bs = [(1.8, [0.61131619810717852, 0.57736760378564378, 0.61131619810717797])]
def rcap(Blist):
    A,b = gen_inequalities(chq.dim)
#    Aeq,beq = gen_equalities(chq.dim)
    Aeq,beq = gen_equalities(chq.dim,Shapley)
    tc = cvx.matrix([ 0.  ,  0.  ,  0.28,  0.31,  0.  ,  0.68,  0.28,  1.  ])
    print A*tc - b
    print Aeq*tc - beq
#    raw_input("Pause")
    At = cvx.spmatrix([],[],[],(A.size[0],1))
    A = cvx.sparse([[A],[At]])
    print A
    Aeqt = cvx.spmatrix([],[],[],(Aeq.size[0],1))
    Aeq = cvx.sparse([[Aeq],[Aeqt]])
#    beq = cvx.matrix(beq,tc='d')
    print Aeq,beq
    
    
#    def objfunc(v,Blist,xk):
#        d = array([ravel(chq.MobiusB(xk[i])-chq.MobiusB(array(Blist[i][1]))) for i in range(len(xk))])
#        return dot(dot(dot(v,d.T),d),v)
#    def fprime(v,Blist,xk):
#        d = array([ravel(chq.MobiusB(xk[i])-chq.MobiusB(array(Blist[i][1]))) for i in range(len(xk))])
#        return 2*dot(dot(d.T,d),v)
    xk = array([1./3,1./3,1./3])
    vk = array([0.,0,0,0,0,0,0,1])
    objfunc = cvx.matrix([0.,0,0,0,0,0,0,0,1])
    print "START"
    R = array([chq.Choquet(array(chq.max_choquet(vk,B[0])),vk) - chq.Choquet(array(B[1]),vk) for B in Blist])
    print "R", R
    print "RSQ", sum(R**2)
    dist = 0
    it = 0
    for i in range(100):
        xk = array(chq.max_choquet(vk,Blist[0][0]))
        diff = chq.Choquet(xk,vk) - chq.Choquet(array(Blist[0][1]),vk)
#        print "XK", xk,array(Blist[0][1])
#        print "DIFF", diff
#        print "CAP", vk
#        print chq.Choquet(xk,vk),chq.Choquet(array(Blist[0][1]),vk)
        if diff < 0:
            print "Bad Optimum!"
            xk = array(Blist[0][1])
            diff = 0
        if diff > dist:
            ZZr = cvx.matrix([chq.MobiusB(xk)-chq.MobiusB(Blist[0][1]), -1]).T
            A = cvx.sparse([A,ZZr])
            b = cvx.matrix([b,0])
#            print "XK", xk
#            print "XR", Blist[0][1]
            sol = cvx.solvers.lp(objfunc,A,b,Aeq,beq,'glpk')
            vk = ravel(sol['x'])[0:8]
            dist = sol['primal objective']
            print "DIST", dist
            print vk
            it = it+1
        else:
            print "FINITA!", it
            break
    print "Capacity", vk
    print "Distance", dist
    print chq.max_choquet(vk,Blist[0][0])
    print Blist[0][1]
        
    print "FINISH"
#    print vk
#    print Blist
#    print chq.max_choquet(vk)
#    R = array([chq.Choquet(array(chq.max_choquet(vk,B[0])),vk) - chq.Choquet(array(B[1]),vk) for B in Blist])
#    Rx = array([array(chq.max_choquet(vk,B[0])) - array(B[1]) for B in Blist])
#    print R
#    print sum(R**2)
#    print Rx
#    print "maxim", chq.Choquet(array(chq.max_choquet(vk,B[0])),vk)
#    print "xr", chq.Choquet(array(B[1]),vk)
#    print sum(array(chq.max_choquet(vk,B[0])))
    
#    print xr 

rcap(Bs)