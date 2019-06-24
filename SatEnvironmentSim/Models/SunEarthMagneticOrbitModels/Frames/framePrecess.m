%Frame Transformation
% function framePrecess
%
% Created: 31.01.2015 10:23:11
% Author: Antoine Pignede
%
% This function calulates the transformation matrix that accounts for the effects
%   of precession.
%
% Copied and edited from Vallado (2013) Fundamentals of Astrodynamics and Applications
%
%  inputs          description                    range / units
%    ttt         - julian centuries of tt
%
%  outputs       :
%    prec        - transformation matrix for mod - j2000 (80 only)
%
%  locals        :
%    convrt      - conversion from arcsec to rad
%    ttt2        - ttt squared
%    ttt3        - ttt cubed
%    zeta        - precession angle               rad
%    z           - precession angle               rad
%    theta       - precession angle               rad
%
%  coupling      :
%    none        -
%
% Note that the original function provides three different formulas (b1950,
%   1980 and 03). The test routines use 1980 which is also the option taken
%   here, the differences are small.
%
% See also
%   frameEcef2eci.m, frameEci2ecef.m  - transformation that need framePrecess
%   frameEci2teme.m, frameTeme2eci.m  - transformation that need framePrecess

function [prec] = framePrecess(ttt) %#codegen

    convrt = pi / (180.0*3600.0);
    ttt2= ttt * ttt;
    ttt3= ttt2 * ttt;

    % find the precession angles
    zeta = 2306.2181*ttt + 0.30188*ttt2 + 0.017998*ttt3;
    theta= 2004.3109*ttt - 0.42665*ttt2 - 0.041833*ttt3;
    z    = 2306.2181*ttt + 1.09468*ttt2 + 0.018203*ttt3;

    zeta = zeta  * convrt; 
    theta= theta * convrt;
    z    = z     * convrt;

    coszeta  = cos(zeta);
    sinzeta  = sin(zeta);
    costheta = cos(theta);
    sintheta = sin(theta);
    cosz     = cos(z);
    sinz     = sin(z);

    % build precession matrix
    prec = eye(3);
    prec(1,1) =  coszeta * costheta * cosz - sinzeta * sinz;
    prec(1,2) =  coszeta * costheta * sinz + sinzeta * cosz;
    prec(1,3) =  coszeta * sintheta;
    prec(2,1) = -sinzeta * costheta * cosz - coszeta * sinz;
    prec(2,2) = -sinzeta * costheta * sinz + coszeta * cosz;
    prec(2,3) = -sinzeta * sintheta;
    prec(3,1) = -sintheta * cosz;
    prec(3,2) = -sintheta * sinz;
    prec(3,3) =  costheta;
end
