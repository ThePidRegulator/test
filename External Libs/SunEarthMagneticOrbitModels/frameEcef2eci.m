%Frame Transformations
% function ecef2eci
%
% Created: 30.01.2015 16:39:02
% Author: Antoine Pignede
%
% This function transforms a vector from the earth fixed (itrf) frame, to
%   the eci mean equator mean equinox (j2000) frame.  the results take into account
%   the effects of precession, nutation, sidereal time, and polar motion.
%
% Copied and edited from Vallado (2013) Fundamentals of Astrodynamics and Applications
%
%  inputs          description                    range / units
%    recef       - position vector earth fixed    km
%    timeUTC     - struct containing:
%      year        - year                           1900 .. 2100
%      mon         - month                          1 .. 12
%      day         - day                            1 .. 28,29,30,31
%      hr          - universal time hour            0 .. 23
%      min         - universal time min             0 .. 59
%      sec         - universal time sec             0.0 .. 59.999
%
%  outputs       :
%    reci        - position vector eci            km
%
%  locals        :
%    deltapsi    - nutation angle                 rad
%    meaneps     - mean obliquity of the ecliptic rad
%    omega       -                                rad
%    prec        - matrix for mod - eci 
%    nut         - matrix for tod - mod 
%    st          - matrix for pef - tod 
%    pm          - matrix for ecef - pef 
%    jdUTC       - julian date                    days from 4713 bc
%    all UTC times are converted to terrestrial time (TT)
%
%  coupling      :
%   framePrecess      - rotation for precession       
%   frameNutation     - rotation for nutation          
%   frameSidereal     - rotation for sidereal time     
%   framePolarm       - rotation for polar motion 
%   timeDatetime2jd   - convert date and time values to julian date
%   timeUtc2tt        - convert UTC to TT
%   timeJd2jc         - convert julian date to julian centuries
%
% See also
%   frameEci2ecef.m  - inverse function

function [reci] = frameEcef2eci(recef,timeUTC) %#codegen

    % find the julian centuries of terrestrial time at given UTC
    jdUTC = timeDatetime2jd(timeUTC);
    timeTT = timeUtc2tt(timeUTC);
    jdTT = timeDatetime2jd(timeTT);
    ttt = timeJd2jc(jdTT);

    % get the transformation matrices
    prec = framePrecess(ttt);
    [deltapsi,meaneps,omega,nut] = frameNutation(ttt);
    st = frameSidereal(jdUTC,deltapsi,meaneps,omega);
    pm = framePolarm(jdUTC);

    % transform ecef to eci
    rpef = pm*recef;
    reci = prec*nut*st*rpef;
end
