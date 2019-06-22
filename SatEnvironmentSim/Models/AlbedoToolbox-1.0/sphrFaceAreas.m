% A = sphereFaceAreaMap(i,j,sy,sx)
% 
% Created: 21.04.2019 16:35:00
% Author: Henrik Rudi Haave
%
% For a sphere divided into fixed degree increments(equirectangular mapping)
% unique surface areas are found. These corresponds to a column of the data array 
% representing such a mapping of a sphere.
% This function is based on Dan D. V. Bhanderi cellarea function 
%
%  Inputs          description                    range / units
%                 
%  Outputs       :
%               A- 
%  Locals        :
%       CONST.EMR-
%              sy-
%              sx-
%            dphi-
%          dtheta-
%            dphi-
%  Coupling      :
%   idx2RadiMap(map)
%  See also      :
%
%
% Notes for the future 
%       package matGeom can probobaly replace this
%
% i,j are the indices of the sy sx REFL matrix.
%


function A = sphrFaceAreas(map,radiMap);

CONST.EMR = 6371.01e3;
[sy,sx] = size(map)

% Grid size
dphi = (pi/sy);
dtheta = (2*pi/sx);

phi = radiMap(1,1:sy);
	
% Diagonal points
phimax = phi + dphi/2;
phimin = phi - dphi/2;

A = CONST.EMR^2*dtheta*(cos(phimin)-cos(phimax));

return
