%Frame Transformations
% function frameNutation
%
% Created: 31.01.2015 10:23:11
% Author: Antoine Pignede
%
% This function calulates the transformation matrix that accounts for the
%   effects of nutation.
%
% Copied and edited from Vallado (2013) Fundamentals of Astrodynamics and Applications
%
%  inputs          description                    range / units
%    ttt         - julian centuries of tt
%
%  outputs       :
%    deltapsi    - nutation angle                 rad
%    meaneps     - mean obliquity of the ecliptic rad
%    omega       - delaunay element               rad
%    nut         - transformation matrix for tod - mod
%
%  locals        :
%    ddpsi       - delta psi correction to gcrf   rad
%    ddeps       - delta eps correction to gcrf   rad
%    iar80       - integers for fk5 1980
%    rar80       - reals for fk5 1980
%    ttt2        - ttt squared
%    ttt3        - ttt cubed
%    l           - delaunay element               rad
%    ll          - delaunay element               rad
%    f           - delaunay element               rad
%    d           - delaunay element               rad
%    deltaeps    - change in obliquity            rad
%
%  coupling      :
%    timeJc2jd        - convert julian century to julian date
%    frameFundarg     - find fundamental arguments
%
% Note that the original function takes the two correction to gcrf factors
%   ddpsi and ddeps as input parameters. The linear interpolation of Kelso
%   and Vallado (2013) Earth Orientation Parameter and Space Weather Data for Flight Operations
%   is used instead.
%
% See also
%   frameEcef2eci.m, frameEci2ecef.m  - transformation that need frameNutation

function [deltapsi,meaneps,omega,nut] = frameNutation(ttt) %#codegen
global DEG2RAD IAR80 RAR80

%% Kelso and and Vallado (2013) Earth Orientation Parameter and Space Weather Data for Flight Operations
    jd = timeJc2jd(ttt);
    MJD = jd - 2400000.5;
    conv = pi / (180*3600); 
    ddpsi = (-8e-6 * MJD + 0.2506) * conv;
    ddeps = (-7e-7 * MJD + 0.022) * conv;

%% Vallado (2013) Fundamentals of Astrodynamics and Applications
    ttt2= ttt*ttt;
    ttt3= ttt2*ttt;
    
    % get the delaunay variables
    [l,l1,f,d,omega] = frameFundarg(ttt);
    
    % find nutation angle and change in obliquity
    deltapsi= 0.0;
    deltaeps= 0.0;
    for i= 106:-1: 1
        tempval= IAR80(i,1)*l + IAR80(i,2)*l1 + IAR80(i,3)*f + ...
                 IAR80(i,4)*d + IAR80(i,5)*omega;
        deltapsi= deltapsi + (RAR80(i,1)+RAR80(i,2)*ttt) * sin( tempval );
        deltaeps= deltaeps + (RAR80(i,3)+RAR80(i,4)*ttt) * cos( tempval );
    end

    % correction to gcrf
    deltapsi = rem( deltapsi + ddpsi, 2.0 * pi );
    deltaeps = rem( deltaeps + ddeps, 2.0 * pi );

    % mean and true obliquity of the ecliptic
    meaneps = -46.8150 *ttt - 0.00059 *ttt2 + 0.001813 *ttt3 + 84381.448;
    meaneps = rem( meaneps/3600.0 ,360.0  );
    meaneps = meaneps * DEG2RAD;    
    trueeps  = meaneps + deltaeps;

    cospsi  = cos(deltapsi);
    sinpsi  = sin(deltapsi);
    coseps  = cos(meaneps);
    sineps  = sin(meaneps);
    costrueeps = cos(trueeps);
    sintrueeps = sin(trueeps);

    % build nutation matrix
    nut = eye(3);
    nut(1,1) =  cospsi;
    nut(1,2) =  costrueeps * sinpsi;
    nut(1,3) =  sintrueeps * sinpsi;
    nut(2,1) = -coseps * sinpsi;
    nut(2,2) =  costrueeps * coseps * cospsi + sintrueeps * sineps;
    nut(2,3) =  sintrueeps * coseps * cospsi - sineps * costrueeps;
    nut(3,1) = -sineps * sinpsi;
    nut(3,2) =  costrueeps * sineps * cospsi - sintrueeps * coseps;
    nut(3,3) =  sintrueeps * sineps * cospsi + costrueeps * coseps;
end
