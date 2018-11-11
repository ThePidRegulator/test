%Time Conversions
% function timeJc2jd
%
% Created: 30.01.2015 14:16:29
% Author: Antoine Pignede
%
% This function finds the julian day given the julian century.
%
% Inputs
%   Julian century of the modified julian date (jd - 2451545.0  )/ 36525;
%
% Outputs
%   Julian date (days from 4713 bc)
%
% See also
%    timeJd2jc.m  - inverse function

function [jd] = timeJc2jd(jc) %#codegen
    jd = jc * 36525.0 + 2451545.0;
end
