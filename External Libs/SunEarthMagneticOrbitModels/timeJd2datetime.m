%Time Conversions
% function timeJd2datetime
%
% Created: 30.01.2015 14:23:30
% Author: Antoine Pignede
%
% This function finds the year, month, day, hour, minute and second
%   given the julian date.
%
% Copied and edited from Vallado (2013) Fundamentals of Astrodynamics and Applications
%
%  inputs          description                    range / units
%    jd          - julian date                    days from 4713 bc
%
%  outputs       :
%    time        - struct containing:
%      year        - year                           1900 .. 2100
%      mon         - month                          1 .. 12
%      day         - day                            1 .. 28,29,30,31
%      hr          - universal time hour            0 .. 23
%      min         - universal time min             0 .. 59
%      sec         - universal time sec             0.0 .. 59.999
%
%  locals        :
%    days        - day of year plus fractional
%                  portion of a day               days
%    tu          - julian centuries from 0 h
%                  jan 0, 1900
%    temp        - temporary real values
%    leapyrs     - number of leap years from 1900
%
%  coupling      :
%    timeDays2datetime  - finds month, day, hour, minute and second given days and year
%
% See also
%   timeDatetime2jd.m  - inverse function

function [time] = timeJd2datetime(jd) %#codegen
    
    % ----------------- find year and days of the year ---------------
    temp   = jd-2415019.5;
    tu     = temp / 365.25;
    time.year   = 1900 + floor( tu );
    leapyrs= floor( ( time.year-1901 )*0.25 );
    %days   = temp - ((time.year-1900)*365.0 + leapyrs ) + 0.00000000001; % nudge by 8.64x10-7 sec to get even outputs
    days   = temp - ((time.year-1900)*365.0 + leapyrs );
    
    % ------------ check for case of beginning of a year -------------
    if days < 1.0
        time.year   = time.year - 1;
        leapyrs= floor( ( time.year-1901 )*0.25 );
        days   = temp - ((time.year-1900)*365.0 + leapyrs );
    end

    % ------------------- find remaining data  -----------------------
    time = timeDays2datetime( time.year,days );
    %time.sec= time.sec - 0.00000086400;
end
