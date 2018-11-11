%Time Conversions
% function timeGstime
%
% Created: 31.01.2015 11:26:48
% Author: Antoine Pignede
%
% This function finds the greenwich sidereal time (iau-82).
%
% Copied and edited from Vallado (2013) Fundamentals of Astrodynamics and Applications
%
%  inputs          description                    range / units
%    jdut1       - julian date of ut1             days from 4713 bc
%
%  outputs       :
%    gst         - greenwich sidereal time        0 to 2pi rad
%
%  locals        :
%    temp        - temporary variable for reals   rad
%    tut1        - julian centuries from the
%                  jan 1, 2000 12 h epoch (ut1)
%
%  coupling      :
%    timeJd2jc   - convert julian date to julian centuries
%
% See also

function gst = timeGstime(jdut1) %#codegen
global TWOPI DEG2RAD

    % julian centuries
    tut1= timeJd2jc(jdut1);

    temp = - 6.2e-6 * tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1  ...
           + (876600.0 * 3600.0 + 8640184.812866) * tut1 + 67310.54841;

    % 360/86400 = 1/240, to deg, to rad
    temp = rem( temp*DEG2RAD/240.0,TWOPI );

    % ------------------------ check quadrants --------------------
    if ( temp < 0.0 )
        temp = temp + TWOPI;
    end

    gst = temp;
end
