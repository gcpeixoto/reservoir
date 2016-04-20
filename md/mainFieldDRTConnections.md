# mainFieldDRTConnections.m

The aim of this script is to compute all the connected components found over the whole SPE Model 2 grid through the DRT values and connectivity criterion. For each DRT, a data structure is stored separately.

## Inputs 

**Connectivity criterion:** 6-neighbour and 26-neighbour criteria are available. However, since the post-processed data are used as input for the CMG simulator, the 26-neighbour criterion is not suitable because the flow is computed face-to-face only, in the light of the Finite Volume Method. With the 26-neighbour criterion, diagonal voxels would be included, but with no effect of flow computation.

**Significative clusters:** the user can choose which is the minimum number of voxels desired per component to save into `.csv` tables. Since a component may be both a large set of voxels (much more than 1) and an isolated voxel (1 voxel only), many components can be computed, thus implying in dozens of files. That is why a limit should be set to quantize what is a significative cluster.

## Strategy

The strategy to compute the DRT network is based on 3 steps:                     

- Find pairs of voxels whose distance of one to another is less than d, where d depends on the connectivity criterion;
- Set up sparse adjacency matrices to create  undirected graphs;
- From the adjacency matrices, all the components are determined by using a connected component algorithm.


