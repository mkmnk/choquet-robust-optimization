'''
Created on 12 Mar 2011

@author: gb
'''


import networkx as nx
from itertools import chain
from time import time
import matplotlib.pyplot as plt
import Choquet_toolpack as chq
from numpy import *

#G = nx.DiGraph()
#G.add_node(1,name='Root',shapley={},int_index={},nec_suff={},func=[],psorted = [11,12,13], cap=array([0, 0.1, 0.1, 0.8, 0.1, 0.8, 0.8, 1]))
#G.add_node(11,name='N1',shapley={},int_index={},nec_suff={},func=[],psorted = [111,112], cap=array([0, 0.8, 0.1, 1]))
#G.add_node(12,name='N2',shapley={},int_index={},nec_suff={},func=[],psorted = [112,133], cap=array([0, 0.3, 0.6, 1]))
#G.add_node(13,name='N3',shapley={},int_index={},nec_suff={},func=[],psorted = [133,134], cap=array([0, 0.7, 0.1, 1]))
#G.add_node(111,name='f1',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x), cap=[])
#G.add_node(112,name='f2',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
#G.add_node(133,name='f3',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
#G.add_node(134,name='f4',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
##
#G.add_edges_from([(11,1),(12,1),(13,1),(111,11),(112,11),(112,12),(133,12),(133,13),(134,13)])
##
##res_vector = {111:0.21, 112:0.25, 133:0.1, 134:0.1}
#f_sorted = [111,112,133,134]


#print chq.max_tree(G,1,f_sorted,1)
#print chq.max_choquet(G.node[11]['cap'],0.35)


#A = time()
#ch_val,weights = node_val(G,1,res_vector)
#ch_weights = [weights[i] for i in f_sorted]
#print ch_val, ch_weights
#print time()-A
#
#
#print G.node



############ MATLAB TEST


#G = nx.DiGraph()
#G.add_node(1,name='Protection',shapley={},int_index={},nec_suff={},func=[],psorted = [2,3], cap=array([0, 0.3, 0.3, 1]))
#G.add_node(2,name='Technical',shapley={},int_index={},nec_suff={},func=[],psorted = [4,5,6], cap=array([0, 0.3, 0.3, 0.6, 0.3, 0.6, 0.6, 1]))
#G.add_node(3,name='Office',shapley={},int_index={},nec_suff={},func=[],psorted = [7,8,9], cap=array([0, 0.2, 0.3, 0.5, 0.3, 0.6, 0.7, 1]))
#G.add_node(4,name='Net. Office',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x), cap=[])
#G.add_node(5,name='FW',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x), cap=[])
#G.add_node(6,name='Part',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
#G.add_node(7,name='NIDS',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
#G.add_node(8,name='IDS',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
#G.add_node(9,name='ENCR',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
##
#G.add_edges_from([(2,1),(3,1),(4,2),(5,2),(6,2),(7,3),(8,3),(9,3)])
##
##res_vector = {111:0.21, 112:0.25, 133:0.1, 134:0.1}
#f_sorted = [4,5,6,7,8,9]



####### AoR EXAMPLE
#
#G = nx.DiGraph()
#G.add_node(1,name='Protection',shapley={},int_index={},nec_suff={},func=[],psorted = [2,21], cap=array([0., 0, 0, 1]))
#G.add_node(2,name='Technical',shapley={},int_index={},nec_suff={},func=[],psorted = [3,17], cap=array([0., 0, 0, 1]))
#G.add_node(3,name='Office',shapley={},int_index={},nec_suff={},func=[],psorted = [4,10,15,16], cap=array([ 0., 0., 0., 0.4, 0., 0., 0.,1.,0.,0.,0.,0.4,0,0.,0.,1.]))
#G.add_node(4,name='Net.Office',shapley={},int_index={},nec_suff={},func=[],psorted = [5,6,7,8,9], cap=array([ 0.  ,  0.29,  0.23,  0.52,  0.  ,  0.33,  0.23,  0.56,  0.  ,0.29,  0.23,  0.52,  0.14,  0.59,  0.37,  0.82,  0.18,  0.47,0.41,  0.7 ,  0.18,  0.51,  0.41,  0.74,  0.18,  0.47,  0.41,0.7 ,  0.32,  0.77,  0.55,  1.  ]))
##
#G.add_node(5,name='FW',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x), cap=[])
#G.add_node(6,name='Part',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
#G.add_node(7,name='NIDS',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
#G.add_node(8,name='IDS',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
#G.add_node(9,name='ENCR',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
##
#G.add_node(10,name='Ser. Office',shapley={},int_index={},nec_suff={},func=[],psorted = [11,12,13,14], cap=array([ 0. ,  0.6,  0. ,  0.6,  0. ,  0.6,  0. ,  0.6,  0. ,  0.6,  0. ,1. ,  0. ,  0.6,  0. ,  1. ]))
##
#G.add_node(11,name='Patch.MGM',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x), cap=[])
#G.add_node(12,name='Log MGM',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x), cap=[])
#G.add_node(13,name='Adm Acc',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
#G.add_node(14,name='Hardening',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
##
#G.add_node(15,name='Phy Offc',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
#G.add_node(16,name='WS Offc',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
##
#G.add_node(17,name='Hosting',shapley={},int_index={},nec_suff={},func=[],psorted = [18,19,20], cap=array([ 0.,  0.,  0.,  0.,  0.,  0.,  1.,  1.]))
##
#G.add_node(18,name='Phy Host',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x), cap=[])
#G.add_node(19,name='Net Host',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
#G.add_node(20,name='Ser Host',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
##
#G.add_node(21,name='Staff',shapley={},int_index={},nec_suff={},func=[],psorted = [22,30], cap=array([0, 1, 0., 1]))
#G.add_node(22,name='Staff. Technical',shapley={},int_index={},nec_suff={},func=[],psorted = [23,29], cap=array([0,1,0., 1]))
#G.add_node(23,name='Access Control',shapley={},int_index={},nec_suff={},func=[],psorted = [24,25,26,27,28], cap=array([ 0.  ,  0.  ,  0.37,  0.6 ,  0.  ,  0.  ,  0.57,  0.8 ,  0.  ,0.  ,  0.37,  0.6 ,  0.  ,  0.  ,  0.77,  1.  ,  0.  ,  0.  ,0.37,  0.6 ,  0.  ,  0.  ,  0.57,  0.8 ,  0.  ,  0.  ,  0.37,0.6 ,  0.  ,  0.  ,  0.77,  1.  ]))
##
#G.add_node(24,name='Sec Policy',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
#G.add_node(25,name='Acc Limit',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
#G.add_node(26,name='User MGM',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x), cap=[])
#G.add_node(27,name='Pol Compl',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x), cap=[])
#G.add_node(28,name='Str auth',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
##
#G.add_node(29,name='Audit',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x), cap=[])
##
#G.add_node(30,name='HR',shapley={},int_index={},nec_suff={},func=[],psorted = [31,32,33], cap=array([ 0.  ,  0.27,  0.5 ,  0.77,  0.23,  0.5 ,  0.73,  1.  ]))
##
#G.add_node(31,name='Sec.Chck',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x), cap=[])
#G.add_node(32,name='Training',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
#G.add_node(33,name='Motivation',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
##
#


###### Thesis EXAMPLE

G = nx.DiGraph()
G.add_node(1,name='Protection',shapley={},int_index={},nec_suff={},func=[],psorted = [2,21], cap=array([0., 0, 0, 1]))
G.add_node(2,name='Technical',shapley={},int_index={},nec_suff={},func=[],psorted = [3,17], cap=array([0., 0, 0, 1]))
G.add_node(3,name='Office',shapley={},int_index={},nec_suff={},func=[],psorted = [4,10,15,16], cap=array([ 0.       , 0.       , 0.       ,  0.4      , 0.       , 0.       , 0.       ,0.4808997, 0.       , 0.       , 0.       ,  0.9191003, 0.       , 0.       ,0.       ,  1.       ]))
G.add_node(4,name='Net.Office',shapley={},int_index={},nec_suff={},func=[],psorted = [5,6,7,8,9], cap=array([ 0.        ,  0.22546259,  0.02067229,  0.45074014, 0.        ,  0.24363598,0.03884569,  0.46891354, 0.        ,  0.22546259,  0.02067229,  0.45074014,0.        ,  0.24363598,  0.03884569,  0.46891354, 0.        ,  0.22546259,0.23316051,  0.66322837,  0.        ,  0.24992341,  0.25762134,  0.68768919,0.        ,  0.39950406,  0.31381892,  0.83726984,  0.        ,  0.48754164,0.33827974,  1.        ]))
#
G.add_node(5,name='FW',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x), cap=[])
G.add_node(6,name='Part',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
G.add_node(7,name='NIDS',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
G.add_node(8,name='IDS',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
G.add_node(9,name='ENCR',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
#
G.add_node(10,name='Ser. Office',shapley={},int_index={},nec_suff={},func=[],psorted = [11,12,13,14], cap=array([ 0.        ,  0.32823425, 0.        ,  0.32823425, 0.        ,  0.32823425,0.        ,  0.32823425, 0.        ,  0.32823425, 0.        ,  0.78238592,0.        ,  0.54584833, 0.        ,  1.        ]))
#
G.add_node(11,name='Patch.MGM',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x), cap=[])
G.add_node(12,name='Log MGM',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x), cap=[])
G.add_node(13,name='Adm Acc',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
G.add_node(14,name='Hardening',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
#
G.add_node(15,name='Phy Offc',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
G.add_node(16,name='WS Offc',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
#
G.add_node(17,name='Hosting',shapley={},int_index={},nec_suff={},func=[],psorted = [18,19,20], cap=array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.7,  1.]))
#
G.add_node(18,name='Phy Host',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x), cap=[])
G.add_node(19,name='Net Host',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
G.add_node(20,name='Ser Host',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
#
G.add_node(21,name='Staff',shapley={},int_index={},nec_suff={},func=[],psorted = [22,30], cap=array([ 0.        ,  0.71709234, 0.        ,  1.        ]))
G.add_node(22,name='Staff. Technical',shapley={},int_index={},nec_suff={},func=[],psorted = [23,29], cap=array([ 0.  ,0.85062131, 0. ,1.]))
G.add_node(23,name='Access Control',shapley={},int_index={},nec_suff={},func=[],psorted = [24,25,26,27,28], cap=array([ 0.        , 0.        ,  0.29924119,  0.48669421, 0.        , 0.        ,0.34882735,  0.67414796, 0.        , 0.        ,  0.29924119,  0.48669421,0.        , 0.        ,  0.40562764,  0.85521453, 0.        , 0.        ,0.29924119,  0.48669421, 0.        , 0.        ,  0.34882735,  0.67414796,0.        , 0.        ,  0.29924119,  0.48669421, 0.        , 0.        ,0.40562764,  1.        ]))
#
G.add_node(24,name='Sec Policy',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
G.add_node(25,name='Access Limit',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
G.add_node(26,name='User MGM',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x), cap=[])
G.add_node(27,name='Pol Compl',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x), cap=[])
G.add_node(28,name='Str auth',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
#
G.add_node(29,name='Audit',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x), cap=[])
#
G.add_node(30,name='HR',shapley={},int_index={},nec_suff={},func=[],psorted = [31,32,33], cap=array([ 0.        , 0.        ,  0.37948891,  0.56025555,  0.19801753,  0.43974445,0.75827307,  1.        ]))
#
G.add_node(31,name='Sec.Chck',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x), cap=[])
G.add_node(32,name='Training',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
G.add_node(33,name='Motivation',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])


#
#print G.node[4]['psorted']
#print [[(p,q) for p in G.node[q]['psorted']] for q in G]
[[G.add_edge(p,q) for p in G.node[q]['psorted']] for q in G]
#
#res_vector = {111:0.21, 112:0.25, 133:0.1, 134:0.1}
f_sorted = [5,6,7,8,9,11,12,13,14,15,16,18,19,20,24,25,26,27,28,29,31,32,33]




#chepoints = [ 1.97492791,  1.78183148,  1.43388374,  1.        ,  0.56611626,  0.21816852]
#
#chevalues = [chq.max_tree(G,22,f_sorted,i) for i in chepoints]
#chevalues.reverse()
#print chevalues

print chq.max_tree(G,1,f_sorted,8)
#print chq.max_choquet_sbgr(G,1,f_sorted)
#print chq.max_choquet(G.node[11]['cap'],0.35)
#print G.edges()
#nx.draw_graphviz(G)
#plt.show()


# NODE 4 APPROxIMATION

f4 = [-0.012172192856151,   0.085769178985166,  -0.240455261700179,   0.358000754930389,  -0.405480456240304,   0.670533309153775, 0]
f10 = [-0.023542563307189,   0.165035790762237,  -0.459447035593820,   0.672739534617488, -0.692365235426778,   0.872108848505457, 0]
f23 = [-0.033594756463189,   0.237604817402191,  -0.665360513230697,   0.952007278026191, -0.834748697858656,   0.811410443646072, 0]
#f22 = [-0.029335929406647,   0.207346406139707,  -0.579333438490725,   0.821527980327791, -0.698397454788673,   0.687228652793793, 0]
f30 = [-0.027728205126421,   0.197460701672059,  -0.567139344258168,   0.888084577320966, -1.002827533837975,   1.152014732068982, 0]

        
