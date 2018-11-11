%Time Conversions
% function timeDatetime2years
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
%    time        - struct containing:
%      year        - year                           1900 .. 2100
%      mon         - month                          1 .. 12
%      day         - day                            1 .. 28,29,30,31
%      hr          - universal time hour            0 .. 23
%      min         - universal time min             0 .. 59
%      sec         - universal time sec             0.0 .. 59.999
%
%  outputs       :
%    years       - fractional year                  1900.0 .. 2100.999
%
%  coupling      :
%    timeDatetime2days  - find day of the year plus fraction
%
% See also
%   timeYears2datetime.m  - inverse function

function [years] = timeDatetime2years(time) %#codegen

    % check if leap year or not
	leap = 0;
	if (rem(time.year, 4) == 0)
		if ((rem(time.year, 100) == 0) && (rem(time.year, 400) ~= 0))
            leap = 0;
        else
            leap = 1;
        end
    end
	% find day of year plus fraction
	days = timeDatetime2days(time);
	% convert to years plus fraction
	years = time.year + (days - 1)/(365.0 + leap); % "(days - 1)" because days begins with 1
end
