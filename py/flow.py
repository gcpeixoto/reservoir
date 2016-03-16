import networkx as nx
import numpy as np

# getting edge list from file to assemble the directed graph
fid = open('../tmp/capacity','rb')
G = nx.read_edgelist(fid, nodetype=int, create_using = nx.DiGraph(), data=(('deltaP',float),))

# getting source and sink nodes from file
aux = np.loadtxt('../tmp/sinksource',dtype=int)
source = aux[0]
sink = aux[1] 

# max flow network problem based on pressure gradient 
max_flow_value, max_flow = nx.maximum_flow(G,source,sink,capacity='deltaP')
print 'Max flow value = %f ' % max_flow_value
