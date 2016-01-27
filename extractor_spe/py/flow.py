import networkx as nx
import matplotlib.pyplot as plt

fid = open('../tmp/capacity','rb');
G = nx.read_edgelist(fid, nodetype=int, create_using = nx.DiGraph(), data=(('deltaP',float),))

aux = G.nodes()
source = aux[0]
sink = aux[20]

max_flow = nx.maximum_flow_value(G,source,sink,capacity='deltaP')
print 'Maximum Flow Value = ' + str( max_flow )

#print G.nodes()
#print G.edges()
#print G.edges(data=True)

#nx.draw(G)
#plt.show()
