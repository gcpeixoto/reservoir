%% mainMRSTThreePhaseBlackOil.m - MRST-based 3-phase black oil simulator
%   author:  Dr. Gustavo Peixoto de Oliveira
%            @Federal University of Paraiba
%   mail: gustavo.oliveira@ci.ufpb.br    
%   date: Sep 12nd, 2016        
%             
%   description: performs black-oil simulations over clusters domains                


clear all, close all
mrstModule add ad-core ad-blackoil ad-props mrst-gui spe10

%% Load SPE10 data

rock = SPE10_rock();
rock.perm = convertFrom(rock.perm, milli*darcy);

%% SPE10 grid
layers = 1:85;
cartDims = [60, 220, numel(layers)];
%{
  Computing real dimensions

  Transforms to original measures (ft) and, finally, converts to m.
  1 ft = 0.3048 m
%}
physDims = cartDims.*[20,10,2]*ft*0.3048;

% checks if G.mat grid of SPE10 exists or not
if exist('../mat/G.mat','file')
    load('../mat/G.mat','G');
else
    G = computeGeometry(cartGrid(cartDims, physDims));    
    save('../mat/G.mat','G');
end


%% Input

drtVal = 9;
nc = 1;

%% Load files 

% metrics data structure
load( strcat('../mat/Field/DRT_',num2str(drtVal),'_MetricsData.mat'),'metrics' );
load( strcat('../mat/Field/DRT_',num2str(drtVal),'_LinRegrData.mat'),'linregr' );

% DRT data structure
load( strcat('../mat/Field/DRT_',num2str(drtVal),'.mat'),'drtSt' );

%% Cluster domain

cvi = drtSt.compVoxelInds{ metrics.idComp{nc} }; 
clo = metrics.closenessCentrality{nc};
maxclo = max(clo);
imaxclo = find(clo == maxclo);  
imaxclo = imaxclo(1);

Gc = extractSubgrid(G,cvi);       

rockC.poro = rock.poro(cvi);
rockC.perm = rock.perm(cvi,:);

% normalized porosity
rockC.kn = sqrt(rockC.perm(:,1).^2 + rockC.perm(:,2).^2 + rockC.perm(:,3).^2);
%plotCellData(Gc,rockC.kn);
%% Define well and simulation schedule

simTime = 20*year;
nstep   = 15;
refine  = 5;

rate = 3000;         % IMEX STL = 3000 [m3/day]
bhp = 1.96133*barsa; % IMEX BHP = 2 [kg/cm2]
radius = 0.0762;     % IMEX radius = 0.0762 [m]
refDepth = 480;      % IMEX REFDEPTH = 480 [m];

% Producers
pv = poreVolume(Gc, rockC);
W = addWell([], Gc, rockC, imaxclo,...
            'Name', strcat('Producer M',num2str(nc)),...            
            'Radius',radius,...
            'Val', bhp,...
            'InnerProduct','ip_tpf',...
            'skin',0.0,...
            'Type', 'bhp',...
            'Dir','z',...            
            'refDepth',refDepth,...
            'Sign',-1);  
        
W = addWell(W, Gc, rockC, 600,...
            'Name', strcat('Producer M',num2str(nc+1)),...            
            'Radius',radius,...
            'Val', bhp,...
            'InnerProduct','ip_tpf',...
            'skin',0.0,...
            'Type', 'bhp',...
            'Dir','z',...            
            'refDepth',refDepth,...
            'Sign',-1);    
        
W = addWell(W, Gc, rockC, 1000,...
            'Name', strcat('Producer M',num2str(nc+2)),...            
            'Radius',radius,...
            'Val', bhp,...
            'InnerProduct','ip_tpf',...
            'skin',0.0,...
            'Type', 'bhp',...
            'Dir','z',...            
            'refDepth',refDepth,...
            'Sign',-1);            
        
        

% Compute the timesteps
startSteps = repmat((simTime/(nstep + 1))/refine, refine, 1);
restSteps =  repmat(simTime/(nstep + 1), nstep, 1);
timesteps = [startSteps; restSteps];

% Set up the schedule containing both the wells and the timesteps
schedule = simpleSchedule(timesteps, 'W', W);

%% Set up simulation model

% Three-phase template model
%{
    props @ref. pressure IMEX REFPRES = 48 [kg/cm2];
                         IMEX REFPW = 50 [kg/cm2];
    
    From Edson's .dat PVT table, interpolating;
        
%}
pkgcm2 = 48;
P = [43.1488,48.4133];
BO = [1.07874,1.08603];
BG = [1/41.4432,1/47.2006];
viso = [3.57539,3.32665];
visg = [0.0122167,0.0124432];

muW  = 0.552917;
muO  = interp1(P,viso,pkgcm2);
muG  = interp1(P,visg,pkgcm2);

bO   = interp1(P,BO,pkgcm2);
bG   = bO; % interp1(P,BG,pkgcm2);
bW   = 1/1.01167; % IMEX BWI = 1.01167

pRef = 47.0719*barsa;

rhoW = 988.71; 
rhoO = 887.0;
rhoG = 0.85;

cR = 0.00196133/barsa;    % IMEX CPOR = 200e-5 [kg/cm2];
cW = 4.301912575/barsa;   % IMEX 4.38673e-5 [kg/cm2];
cO = 0.041844976/barsa;   % IMEX 4.267e-5 [kg/cm2];
cG = cO; % not declared(?) 

fluid = initSimpleADIFluid('phases','WOG',...
                           'mu', [muW,muO,muG]*centi*poise, ...
                           'rho',[rhoW, rhoO, rhoG]*kilogram/meter^3,...
                           'n',  [2, 2, 2],...
                           'pRef', pRef,...
                           'cR',cR);%,...
                           %'b',[bW,bO,bG]);%,...                                                      
                           %'c',[cW,cO,cG]);%,... % addition of this input breaks solution
                           
                           

fluid.bO = @(p) exp((p - pRef)*cO);
%{
 Constant oil compressibility
c        = 0.001/barsa;
p_ref    = 300*barsa;
fluid.bO = @(p) exp((p - p_ref)*c);
clf
p0 = (100:10:500)*barsa;
plot(p0/barsa, fluid.bO(p0))
xlabel('Pressure (bar)')
ylabel('Ratio')
title('Reciprocal formation volume factor for oil (bO)')
%}
                           
% Construct reservoir model
gravity reset on
T  = computeTrans(Gc, rockC);
model = ThreePhaseBlackOilModel(Gc, rockC, fluid);

%% Define initial state

sW = 0.1;
sO = 0.7;
sG = 0.2;

sat = [sW*ones(Gc.cells.num,1),sO*ones(Gc.cells.num,1),sG*ones(Gc.cells.num,1)];

g = model.gravity(3);
% Compute initial pressure
p_res = pRef + g*Gc.cells.centroids(:, 3).*...
   (sW.*model.fluid.rhoWS  + sO*model.fluid.rhoOS + sG.*model.fluid.rhoGS);
state0 = initResSol(Gc, p_res, sat);

%%{ 
%plot
clf
plotCellData(Gc, state0.s(:,1))
plotWell(Gc,W,'height',0.5,'radius',0.5,'color','k')
view(50, 50), axis tight
%}

%% Simulate base case

fn = getPlotAfterStep(state0, model, schedule,'view',[50 50], ...
                     'field','s:1','wells',W);
[wellSols, states, report] = ...
   simulateScheduleAD(state0, model, schedule,'afterStepFn',fn);

[QWS,QOS,QGS,BHP] = wellSolToVector(wellSols);
%% Plot reservoir states
% We launch a plotting tool for the reservoir quantities (pressures
% and saturations, located in states).
% figure(1)
% plotToolbar(Gc, states)
% view(50, 50);
% plotWell(Gc,W);