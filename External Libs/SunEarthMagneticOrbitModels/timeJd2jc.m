%Time Conversions
% function timeJd2jc
%
% Created: 30.01.2015 14:16:29
% Author: Antoine Pignede
%
% This function finds the julian century given the julian date.
%
% Inputs
%   Julian date (days from 4713 bc)
%
% Outputs
%   Julian century of the modified julian date (jd - 2451545.0  )/ 36525;
%
% See also
%    timeJc2jd.m  - inverse function

function [jc] = timeJd2jc(jd) %#codegen
    jc = (jd - 2451545.0) / 36525.0;
end
