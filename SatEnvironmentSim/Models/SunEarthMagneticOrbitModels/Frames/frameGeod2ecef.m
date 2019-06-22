%Frame Transformations
% function frameGeod2ecef
%
% Created: 30.01.2015 15:51:02
% Author: Antoine Pignede
%
% Converts geodetic coordinates LATITUDE, LONGITUDE, and ALTITUDE to
%   Earth-centered, Earth fixed (ECEF) coordinates X, Y, and Z. The World Geodetic
%   System 1984 (WGS84) ellipsoid model of the Earth is assumed.
%
% Copied and edited from Compston (2011) Matlab File Exchange International Geomagnetic Reference Field (IGRF) Model
% 
% Inputs:
%   -rgeod: Geodetic latitude and longitude in degrees, height above the Earth in kilometers.
% 
% Ouputs:
%   -recef: x,y,z coordinates of the point in kilometers.
% 
% See also
%   frameEcef2geod.m  - inverse function

function [recef] = frameGeod2ecef(rgeod) %#codegen
global REQU ECCEARTHSQRD DEG2RAD

    % Split input vector in components
    latitude = rgeod(1);
    longitude = rgeod(2);
    altitude = rgeod(3);

    latitude = latitude* DEG2RAD;
    longitude = longitude* DEG2RAD;

    % Conversion from:
    % en.wikipedia.org/wiki/Geodetic_system#Conversion_calculations
    sinlat = sin(latitude);
    coslat = cos(latitude);
    Nphi = REQU ./ sqrt(1 - ECCEARTHSQRD*sinlat.^2);
    x = (Nphi + altitude).*coslat.*cos(longitude);
    y = (Nphi + altitude).*coslat.*sin(longitude);
    z = (Nphi.*(1 - ECCEARTHSQRD) + altitude).*sinlat;

    % Build output vector
    recef = [x;y;z];
end
