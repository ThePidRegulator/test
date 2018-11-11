%Frame Transformations
% function frameEcef2geod
%
% Created: 30.01.2015 15:51:02
% Author: Antoine Pignede
%
% Convert Earth-centered, Earth fixed (ECEF) coordinates X, Y, and Z to
%   geodetic coordinates LATITUDE, LONGITUDE, and ALTITUDE. The World Geodetic
%   System 1984 (WGS84) ellipsoid model of the Earth is assumed.
%
% Copied and edited from Compston (2011) Matlab File Exchange International Geomagnetic Reference Field (IGRF) Model
% 
% Inputs:
%   -recef: x,y,z coordinates of the point in kilometers.
% 
% Ouputs:
%   -rgeod: Geodetic latitude and longitude in degrees, height above the Earth in kilometers.
%
% Note that two functions are in this file because the latitude needs to
%   be computed recursively, which the recur function does. There is a
%   limit to 20 iterations if the maximum error tolerance in the latitude
%   in radians 1e-12 is not reached (I never had more than 10 iterations).
% 
% See also
%   frameGeod2ecef.m  - inverse function

function [rgeod] = frameEcef2geod(recef) %#codegen
global ECCEARTHSQRD

    % Split input vector in components
    x = recef(1);
    y = recef(2);
    z = recef(3);

    % Longitude is easy
    longitude = atan2(y, x)*180/pi;

    % Compute latitude recursively
    rd = hypot(x, y);
    [latitude, Nphi] = recur(asin(z ./ hypot(x, hypot(y, z))), z, ...
        rd, 1);
    sinlat = sin(latitude);
    coslat = cos(latitude);
    latitude = latitude*180/pi;

    % Get altitude from latitude.
    altitude = rd.*coslat + (z + ECCEARTHSQRD*Nphi.*sinlat).*sinlat - Nphi;

    % Build output vector
    rgeod = [latitude;longitude;altitude];
end

function [latitude, Nphi] = recur(lat_in, z, rd, iter) %#codegen
global REQU ECCEARTHSQRD

    sinlat_in = sin(lat_in);
    thisNphi = REQU ./ sqrt(1 - ECCEARTHSQRD*sinlat_in.^2);
    nextlat = atan((z + thisNphi*ECCEARTHSQRD.*sinlat_in)./rd);
    
    if all(abs(lat_in - nextlat) < 1e-12) || iter > 20
        latitude = nextlat;
        Nphi = thisNphi;
    else
        [latitude, Nphi] = recur(nextlat, z, rd, iter + 1);
    end
end
