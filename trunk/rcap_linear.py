'''
Created on 6 Mar 2011

@author: gb
'''


from sip1 import Shapley,Int_index,Necessity
from capacity_parameters import * 
set_printoptions(edgeitems='Inf',linewidth=90)
cvx.printing.options['dformat'] = '%.3f'
cvx.printing.options['width'] = -1


##### NODE 4
#Bs = [(1, [0.243594872638626, 0.24359487249492975, 0.14555014019949786, 0.14884885968333911, 0.21841125498360731]),
#      (1.3,[0.30327851852498805, 0.30224757973912003, 0.20627495878445778, 0.20627495878459071, 0.28192398416684339]),
#      (1.5,[0.34213054214241889, 0.34206832277528476, 0.23366519105260558, 0.25574834748088104, 0.3263875965488095]),
##      (4, [0.84215360521666827, 0.84215360500890624, 0.71947816570909506, 0.76878785421892948, 0.82742676984640096]),
#      (1.7,[0.38237986113904965, 0.38231750582489638, 0.26351997268166361, 0.30470301456184812, 0.3670796457925421]),
#      (2, [0.44215361805164943, 0.4421536180191214, 0.31947838902990883, 0.36878756441482308, 0.4274268104844971])]

##### NODE 10
#Bs = [(1,[0.35092133635456541, 0.21635955454691852, 0.21635955454900077, 0.21635955454951539]),
#      (1.3, [0.41949386527492133, 0.29504585951640067, 0.29041441569222615, 0.2950458595164519]),
#      (1.5, [0.46717371572648747, 0.34852631435681353, 0.33577365556038802, 0.34852631435631093]),
#      (1.7, [0.51539364332015614, 0.40216451360503075, 0.38027732946919185, 0.40216451360562144]),
#      (2, [0.58780355582994115, 0.4822891961571667, 0.44761805185451409, 0.48228919615837801])]

##### NODE 23
#Bs = [(1, [0.18660250734617084, 0.34365555829315531, 0.18660250709324508, 0.16880119526339932, 0.1143382320040295]),
#      (1.3,[0.25053057965158892, 0.40135393266855413, 0.25053057964757169, 0.23329606538486544, 0.16428884264741997]),
#      (1.5,[0.28939959856792341, 0.44370423042124113, 0.28939959856611336, 0.27822701845752279, 0.19926955398719906]),
#      (1.7,[0.32671634623502255, 0.4876687255476812, 0.32671634623086959, 0.32391954318938204, 0.23497903879704451]),
#      (2, [0.39056362170502223, 0.54477663917113561, 0.39056362167393233, 0.38785847407128893, 0.28623764337862079])]

#### NODE 30
#Bs = [(1, [0.23844277505235181, 0.38111964206696153, 0.38043758288068663]),
#      (1.3,[0.32867542288911922, 0.48596014758920958, 0.48536442952167114]),
#      (1.5,[0.3934142893804537, 0.55356315008857937, 0.55302256053096677]),
#      (1.7,[0.46008095262139193, 0.62022976426166065, 0.61968928311694726]),
#      (2, [0.56008095042841932, 0.72022976535902872, 0.71968928421255218])]

#### NODE 3
#Bs = [(1, [0.46921808703786072, 0.36652042129764673, 0.067253467187705557, 0.097008024476787025]),
#      (1.3,[0.60919644585995125, 0.47878471762088459, 0.087169773886195109, 0.12484906263296915]),
#      (1.5,[0.70229008141215521, 0.55386326328402191, 0.10044492663642898, 0.14340172866739392]),
#      (1.7,[0.79528079201207535, 0.62895610123763535, 0.11377396625398514, 0.16198914049630442]),
#      (2, [0.93473020579507382, 0.74135763472995408, 0.13396190991495402, 0.18995024956001788])]

### Node 22

#Bs = [(1, [0.85361149198708419, 0.14638850801291589]),
#      (1.3,[1.1103246428231022, 0.18967535717689785]),
#      (1.5,[1.2812160542502409, 0.2187839457497591]),
#      (1.7,[1.4526840789224018, 0.24731592107759812]),
#      (2, [1.7094156765802691, 0.29058432341973089])]

### Node 21
Bs = [(1, [0.76891189893524137, 0.23108810106475869]),
      (1.3,[0.99715829542584677, 0.30284170457415333]),
      (1.5,[1.1470693030478198, 0.35293069695218021]),
      (1.7,[1.2967403297082274, 0.40325967029177245]),
      (2, [1.5204948049906577, 0.47950519500934241])]

def rcap(Blist):
    A,b = gen_inequalities(chq.dim,Shapley,Int_index)
    Aeq,beq = gen_equalities(chq.dim,Shapley['Sh_values'],Int_index,Necessity)
#    At = cvx.spmatrix([],[],[],(A.size[0],1))
#    A = cvx.sparse([[A],[At]])
#    Aeqt = cvx.spmatrix([],[],[],(Aeq.size[0],1))
#    Aeq = cvx.sparse([[Aeq],[Aeqt]])
    
    
    objfunc =cvx.matrix(zeros(2**chq.dim)) 
    xk = array([1./chq.dim for i in range(chq.dim)])
    sol = cvx.solvers.lp(objfunc,A,b,Aeq,beq,'glpk')
    vk = ravel(sol['x'])[0:2**chq.dim]
############
    print "START"
    R = array([chq.Choquet(array(chq.max_choquet(vk,B[0])),vk) - chq.Choquet(array(B[1]),vk) for B in Blist])
    print "R", R
    print "RSQ", sum(R**2)
    dist = 0
    it = 0
    objfunc = cvx.matrix([cvx.matrix(zeros(2**chq.dim)),1])
    At = cvx.spmatrix([],[],[],(A.size[0],1))
    A = cvx.sparse([[A],[At]])
    Aeqt = cvx.spmatrix([],[],[],(Aeq.size[0],1))
    Aeq = cvx.sparse([[Aeq],[Aeqt]])
    for i in range(200):
        
        xk1 = array(chq.max_choquet(vk,Blist[0][0]))
        xk13 = array(chq.max_choquet(vk,Blist[1][0]))
        xk15 = array(chq.max_choquet(vk,Blist[2][0]))
        xk17 = array(chq.max_choquet(vk,Blist[3][0]))
        xk2 = array(chq.max_choquet(vk,Blist[4][0]))
        diff = max(chq.Choquet(xk1,vk) - chq.Choquet(array(Blist[0][1]),vk), 
                   chq.Choquet(xk13,vk) - chq.Choquet(array(Blist[1][1]),vk),
                   chq.Choquet(xk15,vk) - chq.Choquet(array(Blist[2][1]),vk),
                   chq.Choquet(xk17,vk) - chq.Choquet(array(Blist[3][1]),vk),
                   chq.Choquet(xk2,vk) - chq.Choquet(array(Blist[4][1]),vk)) 
#        print "XK", xk,array(Blist[0][1])
#        print "DIFF", diff
#        print "CAP", vk
#        print chq.Choquet(xk,vk),chq.Choquet(array(Blist[0][1]),vk)
        if diff < 0:
            print "Bad Optimum!"
            xk1 = array(Blist[0][1])
            xk13 = array(Blist[1][1])
            xk15 = array(Blist[2][1])
            xk17 = array(Blist[3][1])
            xk2 = array(Blist[4][1])
            diff = 0
        if diff > dist:
            print "DIFF", diff
            if chq.Choquet(xk1,vk) - chq.Choquet(array(Blist[0][1]),vk) > dist:
                ZZr1 = cvx.matrix([chq.MobiusB(xk1)-chq.MobiusB(Blist[0][1]), -1]).T
                A = cvx.sparse([A,ZZr1])
                b = cvx.matrix([b,0])
            if chq.Choquet(xk13,vk) - chq.Choquet(array(Blist[1][1]),vk) > dist:
                ZZr13 = cvx.matrix([chq.MobiusB(xk13)-chq.MobiusB(Blist[1][1]), -1]).T
                A = cvx.sparse([A,ZZr13])
                b = cvx.matrix([b,0])
            if chq.Choquet(xk15,vk) - chq.Choquet(array(Blist[2][1]),vk) > dist:
                ZZr15 = cvx.matrix([chq.MobiusB(xk15)-chq.MobiusB(Blist[2][1]), -1]).T
                A = cvx.sparse([A,ZZr15])
                b = cvx.matrix([b,0])    
            if chq.Choquet(xk17,vk) - chq.Choquet(array(Blist[3][1]),vk) > dist:        
                ZZr17 = cvx.matrix([chq.MobiusB(xk17)-chq.MobiusB(Blist[3][1]), -1]).T
                A = cvx.sparse([A,ZZr17])
                b = cvx.matrix([b,0])
            if chq.Choquet(xk2,vk) - chq.Choquet(array(Blist[4][1]),vk) > dist:
                ZZr2 = cvx.matrix([chq.MobiusB(xk2)-chq.MobiusB(Blist[4][1]), -1]).T
                A = cvx.sparse([A,ZZr2])
                b = cvx.matrix([b,0])            
            
#            print "XK", xk
#            print "XR", Blist[0][1]
            sol = cvx.solvers.lp(objfunc,A,b,Aeq,beq,'glpk')
            vk = ravel(sol['x'])[0:2**chq.dim]
            dist = sol['primal objective']
            print "DIST",it, dist
#            print vk
            it = it+1
        else:
            print "FINITA!", it
            break
    print "Capacity", vk
#    print "II2-4", cvx.matrix(int_index(2**chq.dim,(0b00010,0b01000)))*cvx.matrix(vk)
#    print "II2-5", cvx.matrix(int_index(2**chq.dim,(0b00010,0b10000)))*cvx.matrix(vk)
#    print "II1-3", cvx.matrix(int_index(2**chq.dim,(0b00100,0b00001)))*cvx.matrix(vk)
#    print "II2-3", cvx.matrix(int_index(2**chq.dim,(0b00010,0b00100)))*cvx.matrix(vk)
    print "SH1", cvx.matrix(shapley(2**chq.dim,0b00001))*cvx.matrix(vk)
    print "SH2", cvx.matrix(shapley(2**chq.dim,0b00010))*cvx.matrix(vk)
#    print "SH3", cvx.matrix(shapley(2**chq.dim,0b00100))*cvx.matrix(vk)
#    print "SH4", cvx.matrix(shapley(2**chq.dim,0b01000))*cvx.matrix(vk)
#    print "SH5", cvx.matrix(shapley(2**chq.dim,0b10000))*cvx.matrix(vk)
    print "Distance1", chq.Choquet(array(chq.max_choquet(vk,Blist[0][0])),vk) - chq.Choquet(array(Blist[0][1]),vk)
    print "Distance13", chq.Choquet(array(chq.max_choquet(vk,Blist[1][0])),vk) - chq.Choquet(array(Blist[1][1]),vk)
    print "Distance15", chq.Choquet(array(chq.max_choquet(vk,Blist[2][0])),vk) - chq.Choquet(array(Blist[2][1]),vk)
    print "Distance17", chq.Choquet(array(chq.max_choquet(vk,Blist[3][0])),vk) - chq.Choquet(array(Blist[3][1]),vk)
    print "Distance2", chq.Choquet(array(chq.max_choquet(vk,Blist[4][0])),vk) - chq.Choquet(array(Blist[4][1]),vk)
    print "Points1", chq.max_choquet(vk,Blist[0][0]), Blist[0][1]
    print "Points13", chq.max_choquet(vk,Blist[1][0]), Blist[1][1]
    print "Points15", chq.max_choquet(vk,Blist[2][0]), Blist[2][1]
    print "Points17", chq.max_choquet(vk,Blist[3][0]), Blist[3][1]
    print "Points2", chq.max_choquet(vk,Blist[4][0]), Blist[4][1]
        
    print "FINISH"
    return vk
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

#rcap(Bs)

