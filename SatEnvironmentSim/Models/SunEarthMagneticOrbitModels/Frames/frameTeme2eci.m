%Frame Transformations
% function frameTeme2eci
%
% Created: 31.01.2015 15:18:33
% Author: Antoine Pignede
%
% This function transforms a vector from the true equator mean equinox of date,
%   (teme) frame to the mean equator mean equinox (j2000) frame.
%
% Copied and edited from Vallado (2013) Fundamentals of Astrodynamics and Applications
%
%  inputs          description                    range / units
%    rteme       - position vector of date
%                    true equator, mean equinox   km
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
%    prec        - matrix for eci - mod
%    nutteme     - matrix for mod - teme - an approximation for nutation
%    all UTC times are converted to terrestrial time (TT)
%
%  coupling      :
%   framePrecess      - rotation for precession        eci - mod
%   frameTruemean     - rotation for truemean          eci - teme
%   timeDatetime2jd   - convert date and time values to julian date
%   timeUtc2tt        - convert UTC to TT
%   timeJd2jc         - convert julian date to julian centuries
%
% See also
%   frameEci2teme.m  - inverse function

function [reci] = frameTeme2eci(rteme,timeUTC) %#codegen

    % find the julian centuries of terrestrial time at given UTC
    timeTT = timeUtc2tt(timeUTC);
    jdTT = timeDatetime2jd(timeTT);
    ttt = timeJd2jc(jdTT);

    % get the transformation matrices
    prec = framePrecess(ttt);
    nutteme = frameTruemean(ttt);

    % transform teme to eci
    reci = prec * nutteme * rteme;
end
