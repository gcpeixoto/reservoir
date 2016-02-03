import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from mayavi import mlab

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

nd = 3
if nd == 2:
 # drawing 
 pos=nx.spring_layout(G) # positions for all nodes

 # nodes
 nx.draw_networkx(G,pos,arrows=True,with_labels=True,node_size=200)
 plt.show()

elif nd == 3:

 # sample code (still not what I want) 
 # reorder nodes from 0,len(G)-1
 G=nx.convert_node_labels_to_integers(G)
 # 3d spring layout
 pos=nx.spring_layout(G,dim=3)
 # numpy array of x,y,z positions in sorted node order
 xyz=np.array([pos[v] for v in sorted(G)])
 # scalar colors
 scalars=np.array(G.nodes())+5

 mlab.figure(1, bgcolor=(0, 0, 0))
 mlab.clf()

 pts = mlab.points3d(xyz[:,0], xyz[:,1], xyz[:,2],
					 scalars,
					 scale_factor=0.1,
					 scale_mode='none',
					 colormap='Blues',
					 resolution=20)

 pts.mlab_source.dataset.lines = np.array(G.edges())
 tube = mlab.pipeline.tube(pts, tube_radius=0.01)
 mlab.pipeline.surface(tube, color=(0.8, 0.8, 0.8))

 #mlab.savefig('flow-graph-3d.png')
 mlab.show() # interactive window
