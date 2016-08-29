
mrstModule add ad-core ad-blackoil ad-props mrst-gui
clear; close all;

%% 

ic = 45; jc = 68;   % center well coordinates
drtValue = 13;      % DRT to choose 
ncluster = 1;

%% LOAD FILES

[PHI,KX,KY,KZ,~,~,~,~,DRT] = loadMatFiles;

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

%% Define wells and simulation schedule

simTime = 20*year;
nstep   = 25;
refine  = 5;

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

CVCFarU = unique(CVCFar,'rows');
id1 = strmatch(CVCFarU(1,:), cvc);
id2 = strmatch(CVCFarU(3,:), cvc);
id3 = strmatch(CVCFarU(4,:), cvc);
cellIndex = [imaxclo;iminclo;id1;id2;id3];

% Producers
pv = poreVolume(G, rock);

% set up the well
rate = 3000*meter^3/day;
bhp  = 1.96133*barsa;
radius = 0.0762; % m

for i = 1:length(cellIndex)

fprintf('Simulating producer %d.\nWell allocated at cell %d',i,cellIndex(i));

W = [];
W = addWell(W, G, rock, cellIndex(i), ...
        'Name','Producer 1', ...
        'Type','bhp','Val', bhp, ... % lin. system breaking with rate
        'Radius',radius,'Skin',0.0,...
        'InnerProduct','ip_tpf','Sign',-1);
    
% Plot the wells
%{
figure(1); plotCellData(G, rock.poro)
hold on
plotWell(G, W, 'radius',0.5,'height',2)
axis tight
%}

% Compute the timesteps
startSteps = repmat((simTime/(nstep + 1))/refine, refine, 1);
restSteps =  repmat(simTime/(nstep + 1), nstep, 1);
timesteps = [startSteps; restSteps];

% Set up the schedule containing both the wells and the timesteps
schedule = simpleSchedule(timesteps, 'W', W);


%% Set up simulation model

% Three-phase template model
fluid = initSimpleADIFluid('phases', 'WOG',...
                            'mu',    [0.553, 3.0, 0.012]*centi*poise, ...
                            'rho',   [988.71, 887, 0.85*1.2222]*kilogram/meter^3, ...
                            'n',     [2, 2, 2],...
                            'pRef', 49.0333*barsa,...
                            'b',[1.01167,1.08545,1.00001],...
                            'c',[(4.38673e-5/0.980665)/barsa,(4.267e-4/0.980665)/barsa,(4.38673e-5/0.980665)/barsa],...
                            'cR',(200e-5/0.980665)/barsa);

% Constant oil compressibility
%c        = (200e-5/0.980665)/barsa;
%c = 0.0001/barsa;
p_ref    = 49.0333*barsa;
%fluid.bO = @(p) exp((p - p_ref)*c);
%clf
%p0 = (10:10:100)*barsa;
%plot(p0/barsa, fluid.bO(p0))
%xlabel('Pressure (bar)')
%ylabel('Ratio')
%title('Reciprocal formation volume factor for oil (bO)')

% Construct reservoir model
gravity reset on
model = ThreePhaseBlackOilModel(G, rock, fluid);


%% Define initial state

sW = 0.001;
sO = 0.998;
sG = 0.001;

sW = sW*ones(G.cells.num,1);
sO = sO*ones(G.cells.num,1);
sG = sG*ones(G.cells.num,1);

sat = [sW,sO,sG];

g = model.gravity(3);

% Compute initial pressure
p_res = p_ref + g*G.cells.centroids(:, 3).*...
    (sW.*model.fluid.rhoWS + sO.*model.fluid.rhoOS + sG.*model.fluid.rhoGS);
state0 = initResSol(G, p_res, sat);

% TODO: plot with interpolation of values of saturations
%figure; 
%plotCellData(G, state0.s(:,1)); % water saturation
%figure; 
%plotCellData(G, state0.s(:,2)); % oil saturation
%figure; 
%plotCellData(G, state0.s(:,3)); % gas saturation

%% Simulate base case

fn = getPlotAfterStep(state0, model, schedule,'view',[50 50], ...
                     'field','s:1','wells',W);
[wellSols, states, report] = ...
   simulateScheduleAD(state0, model, schedule,'afterStepFn',fn);

fname = strcat('states_DRT_',num2str(drtValue),'_Cluster_',num2str(ncluster),'_Case_',num2str(i),'_CellIndex_',num2str(cellIndex(i)));
save(fname,'states');           

end

%% Plot reservoir states
figure(1)
plotToolbar(G, states)
view(50, 50);
plotWell(G,W);


