%Frame Transformation
% function framePolarm
%
% Created: 31.01.2015 11:37:42
% Author: Antoine Pignede
%
% This function calulates the transformation matrix that accounts for polar
%   motion.
%
% Copied and edited from Vallado (2013) Fundamentals of Astrodynamics and Applications
% Merged with IERS BULLETIN-A 12 February 2015 Vol. XXVIII No. 007 found at
%   http://www.iers.org/IERS/EN/Publications/Bulletins/bulletins.html for
%   the approximation of the polar motion coefficients over a long (years) time range.
%
%  inputs          description                    range / units
%    jd          - julian date                    days from 4713 bc
%
%  outputs       :
%    pm          - transformation matrix for ecef - pef
%
%  locals        :
%    convrt      - conversion from arcsec to rad
%    xp          - polar motion coefficient       rad
%    yp          - polar motion coefficient       rad
%
%  coupling      :
%    none.
%
% Note that the original function provides two different formulas (1980 and
%   2000). The test routines use 1980 which is also the option taken here,
%   the differences are small.
%
% See also
%   frameEcef2eci.m, frameEci2ecef.m  - transformation that need framePolarm

function [pm] = framePolarm(jd) %#codegen
global MJDPOLAR XPOLAR YPOLAR

%% IERS BULLETIN-A 12 February 2015 Vol. XXVIII No. 007 http://www.iers.org/IERS/EN/Publications/Bulletins/bulletins.html
    MJD = jd - 2400000.5;
    A = 2*pi*(MJD-MJDPOLAR)/365.25;
    C = 2*pi*(MJD-MJDPOLAR)/435;
    xp = XPOLAR(1) + XPOLAR(2)*cos(A) + XPOLAR(3)*sin(A) + XPOLAR(4)*cos(C) + XPOLAR(5)*sin(C);  
    yp = YPOLAR(1) + YPOLAR(2)*cos(A) + YPOLAR(3)*sin(A) + YPOLAR(4)*cos(C) + YPOLAR(5)*sin(C);

    % convert arcsec to rad
    conv = pi / (180*3600);
    xp = xp * conv;
    yp = yp * conv;

%% Vallado (2013) Fundamentals of Astrodynamics and Applications
    cosxp = cos(xp);
    sinxp = sin(xp);
    cosyp = cos(yp);
    sinyp = sin(yp);

    % build polar motion matrix
    pm = eye(3);
    pm(1,1) =  cosxp;
    pm(1,2) =  0.0;
    pm(1,3) = -sinxp;
    pm(2,1) =  sinxp * sinyp;
    pm(2,2) =  cosyp;
    pm(2,3) =  cosxp * sinyp;
    pm(3,1) =  sinxp * cosyp;
    pm(3,2) = -sinyp;
    pm(3,3) =  cosxp * cosyp;
end
