# mainExtractorSPE.m

This is the primary and underlying script of the suite. Its execution is mandatory to build the SPE Model 2 as Matlab arrays for posterior use in the other script files. 

## Input data

The SPE Model 2 is rebuilt from the bundle of `.dat` files for porosity and permeability available [here](http://www.spe.org/web/csp/datasets/set02.htm). In turn, porosity (phi) and permeability (kx, ky, kz) fields are stored in 3D arrays of size 60 x 220 x 85.  


## Run options

### Well extraction

Two options for well extraction are available for analysis: 
 
- _N wells_: _N_ wells will be chosen randomly for analysis;  
- _1 specific well_: the user inputs the well's surface coordinates _(x,y)_.

### 3D field plots

An option is available to plot the 3D view of the entire field.

### Dispersion plots

Dispersion plots showing the distribution of permeability and porosity versus depth are available

### VTK exporting

Porosity and permeability fields can be exported for `.vtk` for visualization and data analysis.

### CSV exporting 

Raw data for data analysis can be exported to `.csv` files. For each well, a big table is saved to `../csv` having the following information: 

- **i,j,k**: voxel triplet indices
- **phi_e**: effective porosity
- **kx**: permeability - x 
- **ky**: permeability - y 
- **kz**: permeability - z 
- **kn**: norm-2 permeability
- **phiz**: normalized porosity (PMR: pore-matrix ratio)
- **RQI/FZI**: [^1] 
 - **RQI**: Reservoir Quality Index (per voxel)
 - **FZI**: Flow Zone Indicator (per voxel)  
- **DRT**: Discrete Rock Type (per voxel) [^2]
- **FZI/DRT (IG)**: [^3] 
 - **MFZI**: Modified Flow Zone Indicator (per voxel)
 - **MDRT**: Modified Discrete Rock Type (per voxel)
- **FZI***: [^4]
- **Logs**:
 -  **LogPhiz**: log(phiz)
 -  **LogRQI**: log(RQI)

### MAT saving

`.mat` files are stored. 

### Linear regression fit

Automatized FZI best-fit lines through log(phiz) x log(RQI) plot.

### TODO

* Connate water 
* Cementation factor
* Bland-Altman plots

### References

[^1]: AMAEFULE, J.O. et al. Enhanced reservoir description: using core and log data to identify hydraulic (flow) units and predict permeability in uncored intervals/wells. SPE-26436-MS, 1993.
[^2]: GUO, G. et al. Rock typing as an effective tool for permeability and water-saturation modeling: a case study in a clastic reservoir in the Oriente Basin. SPE 97033-PA, 2005.
[^3]: IZADI, M. et al. New Approach in Permeability and Hydraulic-Flow-Unit Determination, SPE Reservoir Evaluation \& Engineering 16 (3), 257-264, 2013.
[^4]: MIRZAEI-PAIAMAN, A. et al. Improved Method to Identify Hydraulic Flow Units for Reservoir Characterization. Energy Technology, 3, 726-733, 2015.