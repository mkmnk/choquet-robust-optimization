import ConfigParser
import re
import math

def parse_node_config(filename,node):
    """
    Inputs: config filename, section name (e.g. 'node1'). Returns lists of values suitable for passing to capacity_parameters functions
    """
    node = str(node)                      # for integer nodes
    config = ConfigParser.ConfigParser()
    config.read(filename)
    sh_order=[]
    sh_equal=[]
    ii_positive=[]
    ii_negative=[]
    ii_order=[]
    ii_equal=[]
    nes=[]
    node_params={}
    if config.has_option(node,'shapley_order'):
        sh_str = config.get(node,'shapley_order')
        m=re.findall('([0-9]+)<([0-9]+)',sh_str)
        if m:
            for i in m:
                sh_order.append((int(i[0]),int(i[1])))
            node_params['Sh_order']=sh_order    
        m=re.findall('([0-9]+)~([0-9]+)',sh_str)
        if m:
            for i in m:
                sh_equal.append((int(i[0]),int(i[1])))
            node_params['Sh_equal']=sh_equal    
        if config.has_option(node,'shapley_delta'):
            node_params['Sh_delta']=float(config.get(node,'shapley_delta')) 
        else:
            node_params['Sh_delta']=0.05   # TBC2 fix this later to use "global" defaults
    if config.has_option(node,'ii_order'):
        ii_str = config.get(node,'ii_order')
        m=re.findall('\(([0-9]+),([0-9]+)\)>0',ii_str)
        if m:
            for i in m:
                ii_positive.append((int(i[0]),int(i[1])))
            node_params['ii_positive']=ii_positive
        m=re.findall('\(([0-9]+),([0-9]+)\)<0',ii_str)
        if m:
            for i in m:
                ii_negative.append((int(i[0]),int(i[1])))
            node_params['ii_negative']=ii_negative    
        m=re.findall('\(([0-9]+),([0-9]+)\)<\(([0-9]+),([0-9]+)\)',ii_str)
        if m:
            for i in m:
                ii_order.append(((int(i[0]),int(i[1])),(int(i[2]),int(i[3]))))
            node_params['ii_order']=ii_order
        m=re.findall('\(([0-9]+),([0-9]+)\)~\(([0-9]+),([0-9]+)\)',ii_str)
        if m:
            for i in m:
                ii_equal.append(((int(i[0]),int(i[1])),(int(i[2]),int(i[3]))))
            node_params['ii_equal']=ii_equal
        if config.has_option(node,'ii_delta'):
            node_params['ii_delta']=float(config.get(node,'ii_delta')) 
        else:
            node_params['ii_delta']=0.05   # TBC2 fix this later to use "global" defaults
    if config.has_option(node,'necessity'):
        ne_str = config.get(node,'necessity')
        m=re.findall('([0-9]+)',ne_str)
        for i in m:
            nes.append(int(i))
        node_params['necessity']=nes
    if config.has_option(node,'function'):
        f_str = config.get(node,'function')
        node_params['function'] = eval("lambda x:" + f_str, math.__dict__)
    if config.has_option(node,'fgrad'):
        f_str = config.get(node,'fgrad')
        node_params['fgrad'] = eval("lambda x:" + f_str, math.__dict__)    
    if config.has_option(node,'subnodes'):
        node_params['subnodes'] = [int(i) for i in config.get(node,'subnodes').split(',')]
    if config.has_option(node,'utility'):
        ut_str = config.get(node,'utility')
        m=re.findall('([0-9\.]+),([0-9\.]+)',ut_str)
        node_params['limits'] = int(m[0][0])
    return node_params

def config_nodes(filename):
    """
    returns all sections (i.e. nodes) for a specific config file
    """
    config = ConfigParser.ConfigParser()
    config.read(filename)
    return config.sections()

def has_capacity(filename,node):
    """
    checks if a certain node has capacity
    """
    config = ConfigParser.ConfigParser()
    config.read(filename)
    return config.has_option(node,'capacity')
    
def relabel(node_params):
    """
    relabels node's subnodes to 1,2,... so that config could be passed to capacity_paramters routines (gen_inequalities and gen_equalities)
    the way it is done (string replace()) is apparently very dangerous and bug-prone. Probably replace later with something better (regexp) TBC3
    """
    sbnode_map = list(enumerate(sorted(node_params['subnodes']),start=1)) # creates a list of pairs (index,element), which acts as (to,from)
    # sbnode_map_binary = [(bin(i),bin(j)) for i,j in sbnode_map]   # and its binary representation
    for i in node_params:
        if i not in ['function','fgrad','criteria_functions','criteria_fgrad','Sh_delta','ii_delta','limits']:
            for smap in sbnode_map:
                node_params[i] = eval(repr(node_params[i]).replace(str(smap[1]),str(smap[0])))
        elif i in ['criteria_functions','criteria_fgrad','limits']:
            node_params[i]=dict((t,node_params[i][f]) for t,f in sbnode_map)
    return node_params

    


