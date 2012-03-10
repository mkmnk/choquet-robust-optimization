import ConfigParser
import re
import math

def parse_node_config(filename,node):
    """
    Inputs: config filename, section name (e.g. 'node1'). Returns lists of values suitable for passing to capacity_parameters functions
    """
    config = ConfigParser.ConfigParser()
    config.read(filename)
    print config
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
        for i in m:
            sh_order.append((bin(int(i[0])),bin(int(i[1]))))
        node_params['Sh_order']=sh_order    
        m=re.findall('([0-9]+)~([0-9]+)',sh_str)
        for i in m:
            sh_equal.append((bin(int(i[0])),bin(int(i[1]))))
        node_params['Sh_equal']=sh_equal    
    if config.has_option(node,'ii_order'):
        ii_str = config.get(node,'ii_order')
        m=re.findall('\(([0-9]+),([0-9]+)\)>0',ii_str)
        for i in m:
            ii_positive.append((bin(int(i[0])),bin(int(i[1]))))
        node_params['ii_positive']=ii_positive
        m=re.findall('\(([0-9]+),([0-9]+)\)<0',ii_str)
        for i in m:
            ii_negative.append((bin(int(i[0])),bin(int(i[1]))))
        node_params['ii_negative']=ii_negative    
        m=re.findall('\(([0-9]+),([0-9]+)\)<\(([0-9]+),([0-9]+)\)',ii_str)
        for i in m:
            ii_order.append(((bin(int(i[0])),bin(int(i[1]))),(bin(int(i[2])),bin(int(i[3])))))
        node_params['ii_order']=ii_order
        m=re.findall('\(([0-9]+),([0-9]+)\)~\(([0-9]+),([0-9]+)\)',ii_str)
        for i in m:
            ii_equal.append(((bin(int(i[0])),bin(int(i[1]))),(bin(int(i[2])),bin(int(i[3])))))
        node_params['ii_equal']=ii_equal
    if config.has_option(node,'necessity'):
        ne_str = config.get(node,'necessity')
        m=re.findall('([0-9]+)',ne_str)
        for i in m:
            nes.append(int(i))
        node_params['necessity']=nes
    if config.has_option(node,'function'):
        f_str = config.get(node,'function')
        node_params['f'] = eval("lambda x:" + f_str, math.__dict__)
    if config.has_option(node,'subnodes'):
        node_params['subnodes'] = [int(i) for i in config.get('node20','subnodes').split(',')]
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
    



