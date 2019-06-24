%Frame Transformations
% function frameFundarg
%
% Created: 31.01.2015 10:23:11
% Author: Antoine Pignede
%
% This function calulates the delauany variables at given julian centuries.
%
% Copied and edited from Vallado (2013) Fundamentals of Astrodynamics and Applications
%
%  inputs          description                    range / units
%    ttt         - julian centuries of tt
%
%  outputs       :
%    l           - delaunay element               rad
%    ll          - delaunay element               rad
%    f           - delaunay element               rad
%    d           - delaunay element               rad
%    omega       - delaunay element               rad
%
%  locals        :
%    ttt2,ttt3,ttt4  - powers of ttt
%
%  coupling      :
%    none        -
%
% Note that the original function provides four different formulas for each
%   delaunay variable (1980, 1996, 2000b, 2010). The test routines use 1980
%   for ecef2eci/eci2ecef and 2010 (but labelled as 1980) for
%   eci2teme/teme2eci. Here just 2010 is used, the differences are small.
%
% See also
%   frameNutation.m, frameTruemean.m  - function that needs frameFundarg

function [l,l1,f,d,omega] = frameFundarg(ttt) %#codegen
global DEG2RAD

    ttt2= ttt*ttt;
    ttt3= ttt2*ttt;
    ttt4= ttt2*ttt2;

    % ---- determine coefficients for iau 2010 nutation theory ----
    l    =  134.96340251  + ( 1717915923.2178 *ttt + ...
            31.8792 *ttt2 + 0.051635 *ttt3 - 0.00024470 *ttt4 ) / 3600.0;
    l1   =  357.52910918  + (  129596581.0481 *ttt - ...
             0.5532 *ttt2 - 0.000136 *ttt3 - 0.00001149*ttt4 )  / 3600.0;
    f    =   93.27209062  + ( 1739527262.8478 *ttt - ...
            12.7512 *ttt2 + 0.001037 *ttt3 + 0.00000417*ttt4 )  / 3600.0;
    d    =  297.85019547  + ( 1602961601.2090 *ttt - ...
             6.3706 *ttt2 + 0.006593 *ttt3 - 0.00003169*ttt4 )  / 3600.0;
    omega=  125.04455501  + (   -6962890.2665 *ttt + ...
             7.4722 *ttt2 + 0.007702 *ttt3 - 0.00005939*ttt4 )  / 3600.0;

    % ---- convert units to rad
    l    = rem( l,360.0  )     * DEG2RAD;
    l1   = rem( l1,360.0  )    * DEG2RAD;
    f    = rem( f,360.0  )     * DEG2RAD;
    d    = rem( d,360.0  )     * DEG2RAD;
    omega= rem( omega,360.0  ) * DEG2RAD;
end
