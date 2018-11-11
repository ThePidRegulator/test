%Frame Transformations
% function frameSidereal
%
% Created: 31.01.2015 11:17:37
% Author: Antoine Pignede
%
% This function calulates the transformation matrix that accounts for the
%   effects of sidereal time. Notice that deltaspi should not be moded to a
%   positive number because it is multiplied rather than used in a
%   trigonometric argument.
%
% Copied and edited from Vallado (2013) Fundamentals of Astrodynamics and Applications
%
%  inputs          description                    range / units
%    jdutc       - julian centuries of utc        days
%    deltapsi    - nutation angle                 rad
%    meaneps     - mean obliquity of the ecliptic rad
%    omega       - long of asc node of moon       rad
%
%  outputs       :
%    st          - transformation matrix for pef - tod
%
%  locals        :
%    gmst        - mean greenwich sidereal time   0 to 2pi rad
%    ast         - apparent gmst                  0 to 2pi rad
%
%  coupling      :
%    timeGstime      - find gmst
%
% Note that the original function uses ut1 instead of utc to get the gmst.
%   For simplicity ut1 is approximated by utc here. This is legitimiate
%   because abs(ut1 - utc) < 0.9s which is negligible.
%
% See also
%   frameEcef2eci.m, frameEci2ecef.m  - transformation that need frameSidereal

function [st]  = frameSidereal(jdutc,deltapsi,meaneps,omega) %#codegen

    % find mean ast
    gmst= timeGstime( jdutc );
    if (jdutc > 2450449.5 )
        ast= gmst + deltapsi* cos(meaneps) ...
            + 0.00264*pi /(3600*180)*sin(omega) ...
            + 0.000063*pi /(3600*180)*sin(2.0 *omega);
      else
        ast= gmst + deltapsi* cos(meaneps);
    end
    ast = rem (ast,2*pi);
    cosast = cos(ast);
    sinast = sin(ast);

    % build sidereal time matrix
    st = eye(3);
    st(1,1) =  cosast;
    st(1,2) = -sinast;
    st(1,3) =  0.0;
    st(2,1) =  sinast;
    st(2,2) =  cos(ast);
    st(2,3) =  0.0;
    st(3,1) =  0.0;
    st(3,2) =  0.0;
    st(3,3) =  1.0;
end
