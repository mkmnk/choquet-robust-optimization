'''
Created on 17 Jun 2011

@author: gb
'''

import networkx as nx
from itertools import chain
from time import time
import matplotlib.pyplot as plt
import Choquet_toolpack as chq
import operator
from numpy import *
from cap_generator import *


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
G.add_node(21,name='Staff',shapley={},int_index={},nec_suff={},func=[],psorted = [22,31], cap=array([ 0.        ,  0.71709234, 0.        ,  1.        ]))
G.add_node(22,name='Staff. Technical',shapley={},int_index={},nec_suff={},func=[],psorted = [23,30], cap=array([ 0.  ,0.85062131, 0. ,1.]))
#
G.add_node(30,name='Audit',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x), cap=[])
#
G.add_node(31,name='HR',shapley={},int_index={},nec_suff={},func=[],psorted = [32,33,34], cap=array([ 0.        , 0.        ,  0.37948891,  0.56025555,  0.19801753,  0.43974445,0.75827307,  1.        ]))
#
G.add_node(32,name='Sec.Chck',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x), cap=[])
G.add_node(33,name='Training',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
G.add_node(34,name='Motivation',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])


Sols1 = []
Sols8 = []
Sols16 = []
i = 0
Sols20 = []
Varr = cap_generate()
for cap23 in Varr:
    G.add_node(23,name='Access Control',shapley={},int_index={},nec_suff={},func=[],psorted = [24,25,26,27,28,29], cap=array(cap23))
    #
    G.add_node(24,name='Sec Policy',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
    G.add_node(25,name='Access Limit',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
    G.add_node(26,name='User MGM',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x), cap=[])
    G.add_node(27,name='Pol Compl',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x), cap=[])
    G.add_node(28,name='Str auth',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-3*x),gradient=lambda x: 3*exp(-3*x),cap=[])
    G.add_node(29,name='IdM',shapley={},int_index={},nec_suff={},psorted = [],func=lambda x: 1-exp(-x),gradient=lambda x: exp(-x),cap=[])
    [[G.add_edge(p,q) for p in G.node[q]['psorted']] for q in G]
    f_sorted = [5,6,7,8,9,11,12,13,14,15,16,18,19,20,24,25,26,27,28,29,30,32,33,34]
#    print G.edges()
#    nx.draw_graphviz(G)
#    plt.show()
    print i
    Sols1.extend([chq.max_tree(G,1,f_sorted,1)])
    Sols8.extend([chq.max_tree(G,1,f_sorted,8)])
    Sols16.extend([chq.max_tree(G,1,f_sorted,16)])
    Sols20.extend([chq.max_tree(G,1,f_sorted,20)])
#    print Sols8
    i = i+1

print "SOL1"
print Sols1
#print sorted(Sols1,key=operator.itemgetter(1),reverse=True)
print "SOL8"
print Sols8
#print sorted(Sols8,key=operator.itemgetter(1),reverse=True)
print "SOL16"
print Sols16
#print sorted(Sols16,key=operator.itemgetter(1),reverse=True)
print "SOL20"
print Sols20
#print sorted(Sols20,key=operator.itemgetter(1),reverse=True)

#Captest = [0.00, 0.07, 0.08, 0.20, 0.08, 0.20, 0.21, 0.38, 0.08, 0.21, 0.21, 0.38, 0.21, 0.39, 0.39, 0.61, 0.08, 0.21, 0.21, 0.38, 0.21, 0.39, 0.39, 0.61, 0.22, 0.39, 0.39, 0.62, 0.40, 0.62, 0.62, 0.90, 0.50, 0.50, 0.50, 0.55, 0.50, 0.55, 0.55, 0.65, 0.50, 0.55, 0.55, 0.65, 0.55, 0.65, 0.65, 0.80, 0.50, 0.55, 0.55, 0.65, 0.55, 0.65, 0.65, 0.80, 0.55, 0.65, 0.65, 0.80, 0.65, 0.80, 0.80, 1.00]
#[0.00, 0.10, 0.10, 0.24, 0.10, 0.24, 0.25, 0.44, 0.11, 0.25, 0.25, 0.45, 0.25, 0.45, 0.45, 0.70, 0.10, 0.25, 0.25, 0.44, 0.25, 0.45, 0.45, 0.69, 0.26, 0.45, 0.45, 0.70, 0.46, 0.70, 0.70, 1.00, 0.50, 0.50, 0.50, 0.55, 0.50, 0.55, 0.55, 0.65, 0.50, 0.55, 0.55, 0.65, 0.55, 0.65, 0.65, 0.80, 0.50, 0.55, 0.55, 0.65, 0.55, 0.65, 0.65, 0.80, 0.55, 0.65, 0.65, 0.80, 0.65, 0.80, 0.80, 1.00,]