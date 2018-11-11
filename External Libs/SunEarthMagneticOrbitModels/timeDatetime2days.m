%Time Conversions
% function timeDatetime2days
%
% Created: 30.01.2015 14:39:31
% Author: Antoine Pignede
%
% This function finds the fractional days through a year given the year,
%   month, day, hour, minute and second.
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
%    days        - day of year plus fraction of a
%                    day                          days
%
%  locals        :
%    lmonth      - length of months of year
%    i           - index
%
%  coupling      :
%    none.
%
% See also
%   timeDays2datetime.m  - inverse function

function [days] = timeDatetime2days(time) %#codegen

    % --------------- set up array of days in month  --------------
    lmonth = 31*ones(1,12);
    for i= 1 : 12
        if i == 2
            if rem(time.year,4) ~= 0
                lmonth(i)= 28;
            else
                if (rem (time.year,100) == 0) && (rem (time.year,400) ~= 0)
                    lmonth(i)= 28;
                else
                    lmonth(i)= 29;
                end
            end
        elseif i == 4 || i == 6 || i == 9 || i == 11
            lmonth(i)= 30;
        end;
    end

    % find day of the year
    i   = 1;
    days= 0.0;
    while (i < time.mon) && ( i < 12 )
        days= days + lmonth(i);
        i= i + 1;
    end

    % add fractional day
    days= days + time.day + time.hr/24.0 + time.min/1440.0 + time.sec/86400.0;
end
