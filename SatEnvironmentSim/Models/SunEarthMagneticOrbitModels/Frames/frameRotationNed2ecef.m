%Frame Transformations
% function frameRotationNed2ecef
%
% Created: 07.02.2015 10:49:52
% Author: Antoine Pignede
%
% Computes the rotation matrix to transfrom the North East Down (NED) geomagnetic
%   vector into an Earth-centered, Earth fixed (ECEF) geomagnetic vector (X, Y, and Z)
%   given the sin and cos of the latitude and longitude at the location where the
%   geomagnetic vector was calculated.
%
% Copied and edited from Fossen (2011) Handbook of Marine Craft Hydrodynamics and Motion Control
% 
% Inputs:
%   -coslat: cosine of the latitude of the location
%   -sinlat: sine of the latitude of the location
%   -coslon: cosine of the longitude of the location
%   -sinlon: sine of the longitude of the location
% 
% Ouputs:
%   -Ned2ecef: North East Down to ECEF rotation matrix
% 
% See also
%   magIgrf.m, magWmm.m  - functions that use frameRotationNed2ecef

function Ned2ecef =  frameRotationNed2ecef(sinlat, coslat, sinlon, coslon)

    Ned2ecef = [-coslon*sinlat -sinlon -coslon*coslat;-sinlon*sinlat coslon -sinlon*coslat;coslat 0 -sinlat];
end
