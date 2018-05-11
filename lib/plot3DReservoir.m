function plot3DReservoir( I,J,K,field,tit,varargin )
%PLOT3DRESERVOIR plot the reservoir 3D array 

% input:
%       I,J,K: grid bounds
%       field: field to plot (3D array)
%         tit: title
%

if nargin == 4 
    tit = '';
end

% mounting meshgrid for plotting     
i = 1:I;
j = 1:J;
k = K:-1:1; % contrary to show plot upward
[JJ,II,KK] = meshgrid(j,i,k);

disp('Plotting 3D field...')
figure 
h = gca;
set(h,'YDir','reverse'); % z-depth growing downward
slice(JJ,II,KK,field,j,i,k)
title(tit)
        
end

