%Time Conversions
% function timeDatetime2jd
%
% Created: 30.01.2015 14:16:29
% Author: Antoine Pignede
%
% This function finds the julian date given the year, month, day, and time.
%
% Copied and edited from Vallado (2013) Fundamentals of Astrodynamics and Applications
%
%  inputs          description                    range / units
%    time        - struct containing:
%      year        - year                           1900 .. 2100
%      mon         - month                          1 .. 12
%      day         - day                            1 .. 28,29,30,31
%      hr          - universal time hour            0 .. 23
%      min         - universal time min             0 .. 59
%      sec         - universal time sec             0.0 .. 59.999
%
%  outputs       :
%    jd          - julian date                    days from 4713 bc
%
%  locals        :
%    none.
%
%  coupling      :
%    none.
%
% See also
%   timeJd2datetime.m  - inverse function

function jd = timeDatetime2jd(time) %#codegen
        
    jd = 367.0 * time.year  ...
         - floor( (7 * (time.year + floor( (time.mon + 9) / 12.0) ) ) * 0.25 )   ...
         + floor( 275 * time.mon / 9.0 ) ...
         + time.day + 1721013.5  ...
         + ( (time.sec/60.0 + time.min ) / 60.0 + time.hr ) / 24.0;
         %- 0.5 * sign(100.0 * time.year + time.mon - 190002.5) + 0.5;
end
