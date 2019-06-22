%Time Conversions
% function timeDays2datetime
%
% Created: 30.01.2015 14:30:56
% Author: Antoine Pignede
%
% This function converts the day of the year, days, to the equivalent month
%   day, hour, minute and second.
%
% Copied and edited from Vallado (2013) Fundamentals of Astrodynamics and Applications
%
%  inputs          description                    range / units
%    year        - year                           1900 .. 2100
%    days        - julian day of the year         0.0  .. 366.0
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
%    none.
%
% See also
%   timeDatetime2days.m  - inverse function

function [time] = timeDays2datetime(year,days) %#codegen

    % --------------- set up array of days in month  --------------
    lmonth = 31*ones(1,12);
    for i= 1 : 12
        if i == 2
            if rem(year,4) ~= 0
                lmonth(i)= 28;
            else
                if (rem (year,100) == 0) && (rem (year,400) ~= 0)
                    lmonth(i)= 28;
                else
                    lmonth(i)= 29;
                end
            end
        elseif i == 4 || i == 6 || i == 9 || i == 11
            lmonth(i)= 30;
        end;
    end

    dayofyr= floor(days );

    i= 1;
    inttemp= 0;
    while ( dayofyr > inttemp + lmonth(i) ) && ( i < 12 )
        inttemp= inttemp + lmonth(i);
        i= i+1;
    end

    time.year = year;
    % find month and day
    time.mon= i;
    time.day= dayofyr - inttemp;

    % find hours minutes and seconds
    temp= (days - dayofyr )*24.0;
    time.hr  = fix( temp );
    temp= (temp-time.hr) * 60.0;
    time.min = fix( temp );
    time.sec = (temp-time.min) * 60.0;
end
