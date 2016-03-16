function plot3DField( I,J,K,F,tname,varargin )
% PLOT3DFIELD 3D plot of the scalar F(I,J,K) in the reservoir

%   input:
%       I,J,K: 3D grid bounds
%       tname: title (optional)

if nargin == 4
    tname = '';
end

% mounting meshgrid for plotting     
i = 1:I;
j = 1:J;
k = K:-1:1; % contrary to show plot upward
[JJ,II,KK] = meshgrid(j,i,k);

disp('Plotting 3D reservoir...')
figure 
h = gca;
set(h,'YDir','reverse'); % z-depth growing downward
slice(JJ,II,KK,F,j,i,k)
title(tname);

        

end

