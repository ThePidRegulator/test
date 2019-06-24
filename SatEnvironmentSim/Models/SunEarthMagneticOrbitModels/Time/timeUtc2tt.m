%Time Conversions
% function timeUtc2tt
%
% Created: 30.01.2015 17:14:19
% Author: Antoine Pignede
%
% This function converts universal time to terrerstrial time calculating
%   deltaT (deltaT = TT-UT) according to http://eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html
%   assuming the year is between 2005 and 2050.
%
% Inputs
%    timeUTC     - struct containing:
%      year        - year                           1900 .. 2100
%      mon         - month                          1 .. 12
%      day         - day                            1 .. 28,29,30,31
%      hr          - universal time hour            0 .. 23
%      min         - universal time min             0 .. 59
%      sec         - universal time sec             0.0 .. 59.999
%
% Outputs
%    timeTT      - struct containing:
%      year        - year                           1900 .. 2100
%      mon         - month                          1 .. 12
%      day         - day                            1 .. 28,29,30,31
%      hr          - terrerstrial time hour         0 .. 23
%      min         - terrerstrial time min          0 .. 59
%      sec         - terrerstrial time sec          0.0 .. 59.999
%
%  coupling      :
%    timeDays2datetime  - finds month, day, hour, minute and second given days and year
%    timeDatetime2days  - finds days and year given month, day, hour, minute and second
%
% See also
%   timeTt2utc.m  - inverse function

function [timeTT] = timeUtc2tt(timeUTC) %#codegen

    % http://eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html algorithm
    y = timeUTC.year + (timeUTC.mon - 0.5)/12;
    deltaT = 62.92 + 0.32217 * (y-2000) + 0.005589 * (y-2000).^2;
    
    % Days of this year in UTC
    daysUTC = timeDatetime2days(timeUTC);

    % Days of this year in TT
    daysTT = daysUTC + deltaT/(86400);

    % Convert back and adjust if TT is already in the next year
    if (daysTT < 366)
        yearTT = timeUTC.year;
        timeTT = timeDays2datetime(yearTT,daysTT);
    else
        % Check if leap year
        if rem(timeUTC.year,4) ~= 0
            leap = 0;
        else
            if (rem (timeUTC.year,100) == 0) && (rem(timeUTC.year,400) ~= 0)
                leap= 0;
            else
                leap= 1;
            end
        end
        
        if ~leap
            yearTT = timeUTC.year+1;
            timeTT = timeDays2datetime(yearTT,daysTT-365);
        elseif (daysTT < 367)
            yearTT = timeUTC.year;
            timeTT = timeDays2datetime(yearTT,daysTT);
        else
            yearTT = timeUTC.year+1;
            timeTT = timeDays2datetime(yearTT,daysTT-366);
        end
    end
end
