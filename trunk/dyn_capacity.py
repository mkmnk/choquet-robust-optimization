"""
created 06.03.2012
Dynamic robust capacity identification functions
"""

import read_config
import sys
import capacity_parameters
from choquet_nc_opt import convert2kadd,convertnadd,nc_dual
from polyhedron import Vrep, Hrep
import cvxopt as cvx
from numpy import *
#import pdb

def robust_id(config,node):
    """
    Robust identification by solving dual
    inputs: config file name, node number.
    it is assumed that all subnodes have functions
    """
    node_params = read_config.parse_node_config(config,node)
    node_params['criteria_functions'] = {}   # subnode functions called "criteria functions" for passing to capacity_parameters methods
    node_params['criteria_fgrad'] = {}   # subnode gradients called "criteria fgrad" for passing to capacity_parameters methods
    for subnode in node_params['subnodes']: 
        try:
            node_params['criteria_functions'][subnode] = read_config.parse_node_config(config,subnode)['function']
        except:
            print "Not all subnodes have functions! ",subnode
            sys.exit(1)            
        try:
            node_params['criteria_fgrad'][subnode] = read_config.parse_node_config(config,subnode)['fgrad']
        except:
            print "Not all subnodes have gradients! ",subnode
            sys.exit(1)                
    node_params = read_config.relabel(node_params)      # relabel DM information to 1,2,... so that it can be passed to capacity_parameters methods
    A,b = capacity_parameters.gen_inequalities(node_params,convex=0)
    Aeq,beq = capacity_parameters.gen_equalities(node_params,k_additive=2)
    # Convert A,b,Aeq,beq to Ax<=b form, reduce to 2add and calculate vertices
    #
    A = matrix(cvx.matrix([A,Aeq,-Aeq]))
    b = matrix(cvx.matrix([b,beq,-beq]))
    # A = matrix(cvx.matrix([A]))
    # b = matrix(cvx.matrix([b,beq,-beq]))
	# Convert to 2-additive for vertex search simplification
	#
    A, b, bas = convert2kadd(A,b)
    # bas.insert(0,0)                       
    p = Hrep(A, b)
    Vm = array(p.generators)
    print Vm
    print shape(Vm)
    Vm = convertnadd(Vm,bas)
    print Vm
    print shape(Vm)
    nc_dual(Vm,node_params,Budg)

        # cap_list = []
    # for Zi in linspace(0,Z):
    #     cap_list.append((Zi,solve_dual(A,b,Aeq,beq,f,Zi)))
    # vr = linear_nitepolation(cap_list)
    # #
    # # Split vertices into convex and nonconvex and call find_cap for each permutation
    # #


def capacity_id(config,node,id_mode):
    """
    wrapper around different cap_id methods
    also traverses the tree until it is possible to do the identification
    inputs: config file name, node name (e.g. section name), id_mode
    """
    print "capacity_id called"
    for subnode in read_config.parse_node_config(config,node)['subnodes']:
        if not read_config.parse_node_config(config,subnode)['function']:
            if 'subnodes' not in read_config.parse_node_config(config,subnode):
                print "Node ",node," is not an aggregation node and does not have a function defined - check config!"
                sys.exit(1)
            else:
                capacity_id(config,subnode,id_mode)   # recursive call
    if id_mode == 'r':                            # all subnodes have functions
        robust_id(config,node)
        # f_approximate(config,node)        # probably should be called from id function
    elif id_mode == 'e':
        entropy_id(config,node)
        # f_approximate(config,node)        # probably should be called from id function
    elif id_mode == 'v':
        variance_id(config,node)
        # f_approximate(config,node)        # probably should be called from id function
    f_approximate(config,node)        # probably should be called from id function

def main():
    """
    asks for config name, mode, node number
    calls identification functions
    """
    print "Enter config file name"
    config = raw_input("----> ")
    print "Enter mode ([a]ll, [n]ode)"
    mode = raw_input("----> ")
    if mode == 'n':
        print "Enter node number. Available nodes are:"
        print read_config.config_nodes(config)
        node = raw_input("----> ")
        print read_config.parse_node_config(config,node)
        if read_config.has_capacity(config,node):
            print "Node ",node," has capacity"
            print read_config.parse_node_config(config,node)
        elif 'subnodes' not in read_config.parse_node_config(config,node):
            print "Node ",node," is not an aggregation node"
        else:
            while 1:
                print "Select identification mode: [r]obust, max_[e]ntropy, min_[v]ariance"
                id_mode = raw_input("----> ")
                if id_mode in ('r','e','v'):
                    capacity_id(config,node,id_mode)
                    break
                else:
                    print "Input not recognized"
    elif mode == 'a':
        print "Node \t has capacity?:"
        for node in read_config.config_nodes(config):
            print node,"\t",read_config.has_capacity(config,node)
    else:
        print "Input not recognized"
        main()

robust_id('example_node.cfg','20')
        # main()
    
