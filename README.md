# SPEReservoir

A Matlab-based toolbox to study petroleum reservoir characterization based on the [SPE 2 Model](http://www.spe.org/web/csp/datasets/set02.htm) synthetic field.

### Metadata

- Authors: 
 - Gustavo PEIXOTO DE OLIVEIRA, D.Sc.
 - Waldir Leite ROQUE, Ph.D.

## Directory tree

### Base dirs 

* _script_: .m script files liable for the main operations.
* _lib_: .m functions and classes that support the script files.
* _dat_: .dat files provided by the SPE dataset (mandatory). 
* _py_: python files used to interface some operations.
* _cpp_: C++ files used to make an interface with SNAP.  

### Post-processing dirs

* _csv_: .csv files that are produced for data analysis.
* _vtk_: .vtk files for geometry visualization.
* _txt_: .txt files that are produced for specific operations, mainly to
  serve CMG input data.
* _mat_: .mat Matlab files (3D arrays and structures)
* _img_: image files that reproduce reservoir and cluster layers. 
* _figs_: general figures and plots.


### Auxiliary dirs

* _tmp_: temporary files used by interfaces.
* _md_: documentation files.
* _log_: log files stored after script execution.


## Script main files

### Family 1 (reconstruction of the SPE model and volume extraction)

* [mainExtractorSPE](md/mainExtractorSPE.md)
* [mainExtractorSPESubvolume](md/mainExtractorSPESubvolume.md)

### Family 2 (operations over the entire field)

* [mainFieldDRTConnections](md/mainFieldDRTConnections.md)
* [mainFieldDRTGraphMetrics](md/mainFieldDRTGraphMetrics.md)

### Family 3 (operations over a _Volume Of Interest_-VOI within the field)

* [mainVOIDRTConnections](md/mainVOIDRTConnections.md)
* [mainVOIDRTGraphDataWell](md/mainVOIDRTGraphDataWell.md)
* [mainVOIDRTGraphMetrics](md/mainVOIDRTGraphMetrics.md)
* [mainVOIMetricsAnalyzer](md/mainVOIMetricsAnalyzer.md)

### Family 4 (conversion of volumes to images)

* [mainVolume2Image](md/mainVolume2Image.md)
* [mainSubVolume2Image](md/mainSubVolume2Image.md)
* [mainVOICluster2Image](md/mainVOICluster2Image.md)

## Third-party dependencies

* [GIBBON code](http://www.gibboncode.org), by Kevin M. Moerman @MIT
  (REMARK: only the functions `ind2patch` and `logic2subind` are used
  for plotting. However the whole toolbox is recommended.) 
* [SNAP library](http://snap.stanford.edu), by Jure Leskovec @Stanford
  (graph manipulation; see instructions on how to compile it on
  UNIX-based machines)
* [Network components](http://danlarremore.com/), by Daniel Larremore @Harvard (clustering decomposition). Function is also available on Matlab File Exchange [here](http://www.mathworks.com/matlabcentral/fileexchange/42040-find-network-components).
* [MIT-SE Network toolbox](http://strategic.mit.edu/downloads.php?page=matlab_networks), by Buonova G., and de Weck, O.L @MIT Strategic Engineering Research Group (graph manipulation)
* [Networkx](http://networkx.github.io/), by Aric Hagberg, Dan Schult
  and Pieter Swart (REMARK: used in python scripts. THIS IS UNDER
  TESTING!)

# Instructions for compilation and build

- Download the third-party dependencies; 
- Arrange your own directory layout or follow the example below:

	```bash
$USER/Programs/gibbon/       # installation directory of Gibbon
$USER/Programs/matlab-tools/ # ditto Matlab tools
$USER/Programs/snap/         # ditto SNAP
$USER/Programs/networkx/     # ditto Networkx (test; not mandatory)
	```

- Define the environment variable `SNAP_DIR` pointing to SNAP
installation directory, then compile SNAP:_

	```bash 
	echo 'export SNAP_DIR=$USER/Programs/snap' > ~/.bashrc
	cd $SNAP_DIR; make all
	```

- Run the script `configurePaths.m` from inside `/lib` or use the 
Matlab GUI `Set Path` to set the paths to the dependencies.
- Compile `.cpp` files to interface with SNAP:

	```bash
	cd /cpp; make
	```

# Some definitions  

- **voxel:** elementary structure; also called **cell** or **block** 
- **field:** the entire SPE 2 Model, i.e. a grid of 60 x 220 x 85 voxels
- **reservoir:** a subset of the field (VOI) based on
a 3D Moore's neighborhood
- **well:** a subset of a reservoir formed by a column of 85 voxels
- **cluster:** an irregular set of voxels connected by a criterion of neighbourhood selection

# Troubleshooting remarks

- This toolbox is operational on Mac OSX and it is supposed to work on
  UNIX-based platforms. For Windows, some adaptions are still required
  to work and the code was not tested therein.
