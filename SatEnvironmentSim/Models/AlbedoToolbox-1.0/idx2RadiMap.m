% A = idx2RadiMap()
% 
% Created: 21.04.2019 16:35:00
% Author: Henrik Rudi Haave
%
% For a equirectangular map of a sphere the polar coordinates of the 
% centers of bins/elements of the map is found and returned as radiMap.
% This radi map is indexed the same way the equirectangular map is, f.eks:
% radiMap(1,2) gives the polar coordiantes the element map(1,2).
%
% Inputs        :   Description                                 range / units
%            map-  Equirectangular spher mapping, sx > sy    [](syXsx)/
% Outputs       :
%   radiMap(1,:)-  phi, South to respected Z North   [pi-Cy:Cy, zeros]/ radians
%   radiMap(2,:)-  theta, respecting X asis             [-pi+Cx:pi-Cx]/ radians
%               -                                                     /
% Locals        :
%             Cy-  Offset from S..N edges to area centers             / radians
%             Cx-  Offset from W..E edges to area centers             / radians
%             sy-  N(1)..S(sy) map elements/bins/areas                / +integer
%             sx-  W..E map elements/bins/areas                       / +integer
% Coupling      :
% 
% See also      :
%    idx2radi(i,j,sy,sx)
%
% Other notes 
%   package matGeom can probobaly replace this
% i,j are indices coresponding to sy sx.
%

function radiMap = idx2RadiMap(map)
	
        [sy, sx] = size(map);
	radiMap = zeros(2,sx);

        % radi distance from map edges to face centers of the sphere 
        Cx = 2*pi/sx/2; 
        Cy = pi/sy/2;  

        phi = linspace(pi-Cy, Cy, sy);
        theta = linspace(-pi+Cx, pi-Cx, sx );
        radiMap(1,1:sy) = phi;
        radiMap(2,1:sx) = theta;
        

return
