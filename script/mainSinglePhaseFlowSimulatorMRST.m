
mrstModule add incomp
clear; close all;

%% 

ic = 45; jc = 68;   % center well coordinates
drtValue = 13;      % DRT to choose 
ncluster = 4;

%% LOAD FILES

% directory
dbase = strcat( '../mat/Well_I',num2str(ic),'_J',num2str(jc),'/' );

matFiles = dir( strcat(dbase,'VOI_DRT*.mat') ); 
DRTFiles = getDRTFile(matFiles,dbase,drtValue);
load(DRTFiles{1},'VOISt');
load(DRTFiles{2},'metrics');


%% SPE 10 grid 

% Define layers
layers = (1:85).';

% Define grid
cartDims = [  60,  220,   numel(layers)];
physDims = [1200, 2200, 2*cartDims(end)] .* ft();   % ft -> m

rock      = SPE10_rock(layers);
rock.perm = convertFrom(rock.perm, milli*darcy);

G = cartGrid(cartDims,physDims);
G = computeGeometry(G);

%% Cluster geometry and petrophysical properties

in_cells = VOISt.compVoxelInds{ncluster};
out_cells = setdiff(G.cells.indexMap,in_cells);
G = removeCells(G,out_cells); % leaves only cluster cells
%figure, plotGrid(G,'FaceColor','g','EdgeColor','k','EdgeAlpha',0.2); view(3)

% cluster porosity, permeability
cvi = VOISt.compVoxelInds{ncluster};
perm = [rock.perm(cvi,1), rock.perm(cvi,2), rock.perm(cvi,3)];
poro = rock.poro(cvi);
rock = makeRock(G,perm,poro);

%{
figure(1); plotCellData(G, rock.perm(:, 1)/(milli*darcy))
view(50, 50), axis tight
colorbar
title('K_x [mD]')
%}


%% Set fluid data
gravity off
fluid      = initSingleFluid('mu' ,    3*centi*poise     , ...
                             'rho', 887*kilogram/meter^3);
                         
%% Define wells and simulation schedule

% Detecting max/min closeness point
cvc = VOISt.compVoxelCoords{ncluster};
closeness = metrics.closenessCentrality{ncluster};
maxclo = max(closeness);
imaxclo = find(closeness == maxclo);    
cvcmaxclo = cvc(imaxclo,:);

minclo = min(closeness);
iminclo = find(closeness == minclo);    
cvcminclo = cvc(iminclo,:);

% farthest voxels
[D2CFar,CVCFar,ilims,jlims,klims] = getDists2Point(cvc,cvcmaxclo,1);    
compD2CFar = D2CFar;
compCVCFar = CVCFar;
compILims = ilims;
compJLims = jlims;
compKLims = klims;

% set up the well
rate = 10*meter^3/day;
bhp  = 1.96133*barsa; % = 2kg/cm2
radius = 0.0762; % m

W = addWell([], G, rock, imaxclo, ...
        'Name','Producer 1', ...
        'Type', 'bhp' , 'Val', bhp, ...
        'Radius', radius, 'InnerProduct', 'ip_tpf');
    
% Plot the wells
%{
figure(1); plotCellData(G, rock.poro)
hold on
plotWell(G, W, 'radius',0.5,'height',2)
axis tight
%}



%% Compute transmissibilities and solve linear system
% First, we compute one-sided transmissibilities for each local face of
% each cell in the grid. Then, we form a two-point discretization by
% harmonic averaging the one-sided transmissibilities and solve the
% resulting linear system to obtain pressures and fluxes.
T    = computeTrans(G, rock);
rSol = initState(G, W, 47.07192*barsa); % =48kg/cm2
rSol = incompTPFA(rSol, G, T, fluid, 'wells', W);

%%
% Plot the fluxes

figure
cellNo  = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
plotCellData(G, sqrt(accumarray(cellNo,  ...
   abs(convertTo(faceFlux2cellFlux(G, rSol.flux), meter^3/day)))), ...
   'EdgeColor', 'k','EdgeAlpha',0.3);
plotWell(G,W,'height',2,'color','k');
colorbar, axis tight off; view(-80,36)
zoom(1.4);

% figure
% cellNo  = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
% plotCellData(G, sqrt(accumarray(cellNo,  ...
%    abs(convertTo(faceFlux2cellFlux(G, rSol.flux), meter^3/day)))), ...
%    'EdgeColor', 'k','EdgeAlpha',0.1,'FaceAlpha',0.25);
% plotWell(G,W,'height',2,'color','k');
% colorbar, axis tight off; view(-80,36)
% zoom(1.4);

%%
% Plot the pressure distribution
figure
plotCellData(G, convertTo(rSol.pressure(1:G.cells.num), barsa), ...
   'EdgeColor','k','EdgeAlpha',0.1);
plotWell(G, W, 'height', 2, 'color', 'k');
colorbar, axis tight off; view(-80,36)
zoom(1.3);
