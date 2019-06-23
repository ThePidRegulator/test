% GRIDANGLE Calculate angle between two grid index pairs.
%
% rho = gridangle(i1,j1,i2,j2,sy,sx)
%
% $Id: gridangle.m,v 1.3 2006/05/17 14:39:17 danji Exp $

function rho = gridangle(i1,j1,i2,j2,radiMap);

phi1 = radiMap(1,i1);
phi2 = radiMap(1,i2);
theta1 = radiMap(2,j1);
theta2 = radiMap(2,j2);

rho = acos(sin(phi1)*sin(phi2)*cos(theta1-theta2)+cos(phi1)*cos(phi2));

return
