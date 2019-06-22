%Time Conversions
% function timeTt2utc
%
% Created: 30.01.2015 17:14:19
% Author: Antoine Pignede
%
% This function converts universal time to terrerstrial time calculating
%   deltaT (deltaT = TT-UT) according to http://eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html
%   assuming the year is between 2005 and 2050.
%
% Inputs
%    timeTT      - struct containing:
%      year        - year                           1900 .. 2100
%      mon         - month                          1 .. 12
%      day         - day                            1 .. 28,29,30,31
%      hr          - terrerstrial time hour         0 .. 23
%      min         - terrerstrial time min          0 .. 59
%      sec         - terrerstrial time sec          0.0 .. 59.999
%
% Outputs
%    timeUTC     - struct containing:
%      year        - year                           1900 .. 2100
%      mon         - month                          1 .. 12
%      day         - day                            1 .. 28,29,30,31
%      hr          - universal time hour            0 .. 23
%      min         - universal time min             0 .. 59
%      sec         - universal time sec             0.0 .. 59.999
%
%  coupling      :
%    timeDays2datetime  - finds month, day, hour, minute and second given days and year
%    timeDatetime2days  - finds days and year given month, day, hour, minute and second
%
% See also
%   timeUtc2tt.m  - inverse function

function [timeUTC] = timeTt2utc(timeTT) %#codegen

    % http://eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html algorithm
    y = timeTT.year + (timeTT.mon - 0.5)/12;
    deltaT = 62.92 + 0.32217 * (y-2000) + 0.005589 * (y-2000)^2;
    
    % Days of this year in TT
    daysTT = timeDatetime2days(timeTT);

    % Days of this year in UTC
    daysUTC = daysTT - deltaT/(86400);

    % Convert back and adjust if UTC is still in the past year
    if (daysUTC >= 1)
        yearUTC = timeTT.year;
        timeUTC = timeDays2datetime(yearUTC,daysUTC);
    else
        % Check if leap year
        if rem(timeTT.year-1,4) ~= 0
            leap = 0;
        else
            if (rem (timeTT.year-1,100) == 0) && (rem(timeTT.year-1,400) ~= 0)
                leap= 0;
            else
                leap= 1;
            end
        end

        if ~leap
            yearUTC = timeTT.year-1;
            timeUTC = timeDays2datetime(yearUTC,daysUTC+365);
        else
            yearUTC = timeTT.year-1;
            timeUTC = timeDays2datetime(yearUTC,daysUTC+366);
        end
    end
end
