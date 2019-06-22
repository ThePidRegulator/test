%Time Conversions
% function timeYears2datetime
%
% Created: 30.01.2015 14:30:56
% Author: Antoine Pignede
%
% This function converts the fractional year, to the equivalent month
%   day, hour, minute and second.
%
% Inspired by Compston (2011) Matlab File Exchange International Geomagnetic Reference Field (IGRF) Model
%
%  inputs          description                    range / units
%    years       - fractional year                  1900.0 .. 2100.999
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
%    dayofyr     - day of year
%    temp        - temporary extended values
%    inttemp     - temporary integer value
%    i           - index
%    lmonth(12)  - integer array containing the number of days per month
%
%  coupling      :
%    timeDatetime2jd  - find julian date
%    timeDays2datetime  - convert fractional days to mon, day, hr, min, sec
%
% See also
%   timeDatetime2years.m  - inverse function

function [time] = timeYears2datetime(years) %#codegen

    time.year = floor(years);
    leap = 0;
	if (rem(time.year, 4) == 0)
        if ((rem(time.year, 100) == 0) && (rem(time.year, 400) ~= 0))
            leap = 0;
        else
            leap = 1;
        end
    end
    % Convert fractional years to fractional day of the year
    days = 1 + (years - time.year) * (365 + leap);
    time = timeDays2datetime(time.year,days);
end
