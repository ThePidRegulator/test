%Frame Transformations
% function frameEci2teme
%
% Created: 31.01.2015 14:49:06
% Author: Antoine Pignede
%
% This function transforms a vector from the mean equator mean equinox frame
%   (j2000) to the true equator mean equinox of date (teme) frame.
%
% Copied and edited from Vallado (2013) Fundamentals of Astrodynamics and Applications
%
%  inputs          description                    range / units
%    reci        - position vector eci            km
%    timeUTC     - struct containing:
%      year        - year                           1900 .. 2100
%      mon         - month                          1 .. 12
%      day         - day                            1 .. 28,29,30,31
%      hr          - universal time hour            0 .. 23
%      min         - universal time min             0 .. 59
%      sec         - universal time sec             0.0 .. 59.999
%
%  outputs       :
%    rteme       - position vector of date
%                    true equator, mean equinox   km
%
%  locals        :
%    prec        - matrix for eci - mod
%    nutteme     - matrix for mod - teme
%    all UTC times are converted to terrestrial time (TT)
%
%  coupling      :
%   framePrecess      - rotation for precession        eci - mod
%   frameTruemean     - rotation for truemean          mod - teme
%   timeDatetime2jd   - convert date and time values to julian date
%   timeUtc2tt        - convert UTC to TT
%   timeJd2jc         - convert julian date to julian centuries
%
% See also
%   frameTeme2eci.m  - inverse function

function [rteme] = frameEci2teme(reci,timeUTC) %#codegen

    % find the julian centuries of terrestrial time at given UTC
    timeTT = timeUtc2tt(timeUTC);
    jdTT = timeDatetime2jd(timeTT);
    ttt = timeJd2jc(jdTT);

    % get the transformation matrices
    prec = framePrecess(ttt);
    nutteme = frameTruemean(ttt);

    % transform eci to teme
    rteme = nutteme'*prec'*reci;
end
