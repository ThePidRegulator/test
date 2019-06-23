%Sun Vector Prediction
% function sun
%
% Created: 31.01.2015 15:33:09
% Author: Antoine Pignede
%
% This function calculates the normed eci mean equator mean equinox (j2000) position vector
%   of the sun given the julian date.  this is the low precision formula and
%   is valid for years from 1950 to 2050.  accuaracy of apparent coordinates
%   is 0.01  degrees.  notice many of the calculations are performed in
%   degrees, and are not changed until later.  this is due to the fact that
%   the almanac uses degrees exclusively in their formulations.
%
% Copied and edited from Vallado (2013) Fundamentals of Astrodynamics and Applications
%
%  inputs          description                    range / units
%    timeUTC     - struct containing:
%      year        - year                           1900 .. 2100
%      mon         - month                          1 .. 12
%      day         - day                            1 .. 28,29,30,31
%      hr          - universal time hour            0 .. 23
%      min         - universal time min             0 .. 59
%      sec         - universal time sec             0.0 .. 59.999
%
%  outputs       :
%    rsun        - eci position vector of the sun normed to magnitude one
%
%  locals        :
%    meanlong    - mean longitude
%    meananomaly - mean anomaly
%    eclplong    - ecliptic longitude
%    obliquity   - mean obliquity of the ecliptic
%    tut1        - julian centuries of ut1 from
%                  jan 1, 2000 12h
%    ttdb        - julian centuries of tdb from
%                  jan 1, 2000 12h
%    all UTC times are converted to terrestrial time (TT)
%
%  coupling      :
%   timeDatetime2jd   - convert date and time values to julian date
%   timeUtc2tt        - convert UTC to TT
%   timeJd2jc         - convert julian date to julian centuries
%   framePrecess      - rotation for precession
%
% Note that the original function uses ut1 instead of utc to get the gmst.
%   For simplicity ut1 is approximated by utc here. This is legitimiate
%   because abs(ut1 - utc) < 0.9s which is negligible.
% The same holds also for the barycentric dynamical time tdb which is
%   approximated with terrestrial time tt. http://en.wikipedia.org/wiki/Barycentric_Dynamical_Time
%   says that abs(tdb - tt) < 0.002s which is even more negligible.
%   Note that the original function approximates tdb with ut1 which results
%   in larger error: tt - ut = 69s in January 2015.
%
% See also

function [rsun] = sun(timeUTC) %#codegen
global TWOPI DEG2RAD

    % find the julian centuries of terrestrial time at given UTC
    jdUTC = timeDatetime2jd(timeUTC);
    tut1= timeJd2jc(jdUTC);
    timeTT = timeUtc2tt(timeUTC);
    jdTT = timeDatetime2jd(timeTT);
    ttt = timeJd2jc(jdTT);

    % find mean longitude, mean anomaly, ecliptic longitude and mean obliquity
    meanlong= 280.460  + 36000.77*tut1;
    meanlong= rem( meanlong,360.0  );  %deg

    ttdb= ttt;
    meananomaly= 357.5277233  + 35999.05034 *ttdb;
    meananomaly= rem( meananomaly*DEG2RAD,TWOPI );  %rad
    if ( meananomaly < 0.0  )
        meananomaly= TWOPI + meananomaly;
    end

    eclplong= meanlong + 1.914666471 *sin(meananomaly) ...
                + 0.019994643 *sin(2.0 *meananomaly); %deg
    eclplong= rem( eclplong,360.0  );  %deg
    eclplong = eclplong *DEG2RAD;

    obliquity= 23.439291  - 0.0130042 *ttdb;  %deg   
    obliquity= obliquity *DEG2RAD;

    % --------- find magnitude of sun vector, )   components ------
    magr= 1.000140612  - 0.016708617 *cos( meananomaly ) ...
                          - 0.000139589 *cos( 2.0 *meananomaly );    % in au's

    % sun position in MOD
    rsun = [0;0;];
    rsun(1)= magr*cos( eclplong );
    rsun(2)= magr*cos(obliquity)*sin(eclplong);
    rsun(3)= magr*sin(obliquity)*sin(eclplong);

    % transform to eci
    prec = framePrecess(ttt);
    rsun = prec*rsun;
    
    % norm output vector
    rsun = rsun/norm(rsun);
end
