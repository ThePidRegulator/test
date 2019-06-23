% EARTHFOV Field of view on earth by spherical coordinates.
%
% result = earthfov(satsph, refl [, type])
%
% satsph is the vector to the satellite in ECEF frame and spherical
% coordinates. refl is the reflectivity data used for plotting the fov.
% type is 'p' for 3D plot and 's' for spherical. If refl or type are
% unspecified no plot is generated. If refl is not specified, refl must be 
% a vector of the latitudal (refl(1)) and longitudal (refl(2)) resolution 
% of the output data.
%
% $Id: earthfov.m,v 1.15 2006/02/23 08:31:33 danji Exp $
% $Id: earthfov.m,       2019/05/30 00:00:00 Henrik Rudi Haave $

function result = earthfov(satsph,refl,type);

CONST.EMR = 6371.01e3;
CONST.d2r = pi/180;

% LEO shortcut
if satsph(3) < CONST.EMR
	satsph(3) = satsph(3) + CONST.EMR;
end

% Circle value
OUTVALUE = 1;

% Discretization parameters
[sy sx] = size( refl.data );

dx = 2*pi/sx;
dy = pi/sy;
result = zeros( sy, sx );

% Small Circle Center
theta0 = satsph(1);
sinPhi0 = sin( satsph(2) );
cosPhi0 = cos( satsph(2) );

% FOV on earth
%sy phi
%sx theta

rho = acos( CONST.EMR/satsph(3) );

% the speed of this may be further increased by using the properti that (rsat-gridpos)*gridposNormal=0;
% at the satelites fov horizon instead of trig or in other words there is also cone circle mapping that gives fov.
result = acos( (sinPhi0*sin( refl.radiMap(1,1:sy) ))'...
                * cos( theta0-refl.radiMap(2, 1:sx) )...
                + (cosPhi0*cos( refl.radiMap(1, 1:sy) ))' )<=rho;


if nargin > 2 && isfield(refl,'data')
    plot_refl(mask(refl.data,result),type,'no colorbar');
end

return
