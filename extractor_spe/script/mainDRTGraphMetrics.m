%% mainDRTGraphMetrics
%   authors: Dr. Gustavo Peixoto de Oliveira
%            Dr. Waldir Leite Roque
%            @Federal University of Paraiba
%   mail: gustavo.oliveira@ci.ufpb.br    
%   date: Oct 27st, 2015        
%             
%   description: compute metrics for the reservoir networks based on 
%                the DRT field value.
%
%   requirements:
%        - pre-computed .mat files
%        - Matlab third-party additional functions
%

%% DEFAULTS
clear all; close all; clc; format long;
delete('../log/DRTgraphMetrics.log');
diary('../log/DRTgraphMetrics.log');
diary on
%profile on 

aux = load('../mat/DRT_Field.mat');
DRT = aux.DRT;

matFiles = dir('../mat/DRT_4.mat'); 
numfiles = length(matFiles);
metrics = cell(2, numfiles);

for k = 1:numfiles 
  st = load( strcat('../mat/',matFiles(k).name) ); 
  Madj = st.drtSt.allAdjMatrix;
  %betw = node_betweenness_faster(st.drtSt.allAdjMatrix);    
  [deg,~,~] = degrees(Madj);
  metrics{1,k} = deg;
  %C = closeness(Madj);
  
  %metrics{2,k} = C;
    
end

%profsave(profile('info'),'../log/profile_report')
close all
diary off