%Geomagnetic Vector Prediction
% function magIgrf
%
% Created: 06.02.2015 08:26:22
% Author: Antoine Pignede
% 
% Calculates the components of the Earth's magnetic field using the
%   International Geomagnetic Reference Field (IGRF) model.
%   Note that the other parameters the IGRF gives can be computed from BX,
%   BY, and BZ as:
%     Horizonal intensity: hypot(BX, BY) (i.e., sqrt(BX.^2 + BY.^2) )
%     Total intensity: hypot(BX, hypot(BY, BZ))
%     Declination: atan2(BY, BX)
%     Inclination: atan(BZ./hypot(BX, BY))
%
% Copied and edited from Compston (2011) Matlab File Exchange International Geomagnetic Reference Field (IGRF) Model
% 
% Inputs:
%    recef         - position vector eci            km
%    timeUTC       - struct containing:
%      year        - year                           1900 .. 2100
%      mon         - month                          1 .. 12
%      day         - day                            1 .. 28,29,30,31
%      hr          - universal time hour            0 .. 23
%      min         - universal time min             0 .. 59
%      sec         - universal time sec             0.0 .. 59.999
% 
% Outputs:
%    Becef        - geomagnetic field vector in nanotesla (nT) in eci frame
%
% Coupling
%    timeDatetime2years     - convert date and time values to fractional years
%    frameEci2ecef     - transform eci position to ecef position
%    frameEcef2geod    - transform ecef position to geodetic position
%
% Note that the original getigrfcoefs, loadigrfcoefs, frame and time
%    transformation functions have been merged into this igrf function.
%    The original vector input function to calculate the geomagnetic field
%    for more than just one position has been removed.
% 
% See also
%    magWmm.m  - alternative geomagnetic vector predictor

function [Becef] = magIgrf(recef,timeUTC) %#codegen
global IGRFGH IGRFGHSV RMEAN REQU RPOL DEG2RAD

%% Edited excerpt of loadigrfcoefs.m from Compston (2011) Matlab File Exchange International Geomagnetic Reference Field (IGRF) Model

    % Convert time to fractional years.
    time = timeDatetime2years(timeUTC);    
    
    % Interpolate g and h linearly using the values at epoch and secular variation coefficients
    gh = IGRFGH + IGRFGHSV*(time - 2015);
    
    % Highest order
    nmax = sqrt(numel(gh) + 1) - 1;

%% Edited excerpt of igrf.m from Compston (2011) Matlab File Exchange International Geomagnetic Reference Field (IGRF) Model

    % Transform eci input vector to geodetic coordinates
    rgeod = frameEcef2geod(recef);
    latitude = rgeod(1);
    longitude = rgeod(2);
    altitude = rgeod(3);
    
    %%% SPHERICAL COORDINATE CONVERSION %%%
    % Convert from geodetic coordinates latitude, longitude, and altitude to
    % geocentric coordinates r (radius), theta (inclination angle from +z axis),
    % and phi (azimuth angle from +x axis).
    % This method was adapted from igrf11syn, which was a conversion of the
    % IGRF subroutine written in FORTRAN.
    coslat = cos(latitude*DEG2RAD);
    sinlat = sin(latitude*DEG2RAD);
    costheta = cos((90 - latitude)*DEG2RAD);
    sintheta = sin((90 - latitude)*DEG2RAD);
    rho = hypot(REQU*sintheta, RPOL*costheta);
    r = sqrt( altitude.^2 + 2*altitude.*rho + ...
        (REQU^4*sintheta.^2 + RPOL^4*costheta.^2) ./ rho.^2 );
    cd = (altitude + rho) ./ r;
    sd = (REQU^2 - RPOL^2) ./ rho .* costheta.*sintheta./r;
    oldcos = costheta;
    costheta = costheta.*cd - sintheta.*sd;
    sintheta = sintheta.*cd + oldcos.*sd;

    % Convert longitude to radians
    phi = longitude*DEG2RAD;

    % We need cos(m*phi) and sin(m*phi) multiple times
    cosphi = cos((1:nmax)*phi);
    sinphi = sin((1:nmax)*phi);

    %%% BEGIN MAGNETIC FIELD CALCULATION %%%
    % Initialize variables used in for loop below
    Pmax = (nmax+1)*(nmax+2)/2;
    Br = 0; 
    Bt = 0;
    Bp = 0;
    P = zeros(1, Pmax);
    P(1) = 1;
    P(3) = sintheta;
    dP = zeros(1, Pmax); 
    dP(1) = 0; 
    dP(3) = costheta;
    m = 1; 
    n = 0; 
    coefindex = 1;

    % Relative radius of position with respect to mean earth radius squared
    a_r = (RMEAN/r)^2;

    % Increment through all the n's and m's. gh will be a vector with g
    % followed by h for incrementing through n and m except when h would be
    % redundant (i.e., when m = 0).

    for Pindex = 2:Pmax

        % Increment to the next n when m becomes larger than n.
        if n < m
            m = 0;
            n = n + 1;
            a_r = a_r*(RMEAN/r); % We need (Rearth_km./r)^(n+2)
        end

        % Calculate P and dP. They are given recursively according to:
        % 
        % P(0, 0) = 1, P(1, 1) = sin(theta) <- Specified above
        % P(n, n) = sqrt(1 - 1/(2n))*sin(theta)*P(n-1, n-1)
        % P(n, m) = (2n - 1)/sqrt(n^2 - m^2)*cos(theta)*P(n-1, m) -
        %     sqrt(((n-1)^2 - m^2) / (n^2 - m^2)) * P(n-2, m)
        % 
        % dP(0, 0) = 0, dP(1, 1) = cos(theta) <- Specified above
        % dP(n, n) = sqrt(1 - 1/(2n))*(sin(theta)*dP(n-1, n-1) +
        %     cos(theta)*P(n-1, n-1))
        % dP(n, m) = (2n - 1)/sqrt(n^2 - m^2)*(cos(theta)*dP(n-1, m) -
        %     sin(theta)*P(n-1, m)) - sqrt(((n-1)^2 - m^2)/(n^2 - m^2))*
        %     dP(n-2, m)
        if m < n && Pindex ~= 3 % (Pindex=3 is n=1, m=1, initial cond. above)
            last1n = Pindex - n;
            last2n = Pindex - 2*n + 1;
            P(Pindex) = (2*n - 1)/sqrt(n^2 - m^2)*costheta*P(last1n) - ...
                sqrt(((n-1)^2 - m^2) / (n^2 - m^2)) * P(last2n);
            dP(Pindex) = (2*n - 1)/sqrt(n^2 - m^2)*(costheta*dP(last1n) - ...
                sintheta*P(last1n)) - sqrt(((n-1)^2 - m^2) / (n^2 - m^2)) * ...
                dP(last2n);
        elseif Pindex ~= 3
            lastn = Pindex - n - 1;
            P(Pindex) = sqrt(1 - 1/(2*m))*sintheta*P(lastn);
            dP(Pindex) = sqrt(1 - 1/(2*m))*(sintheta*dP(lastn) + ...
                costheta*P(lastn));
        end

        % Calculate the magnetic field components as a running sum. Find
        % explicit expressions for these in Global Earth Physics: a Handbook of
        % Physical Constants by Thomas J. Aherns (1995), pg. 49. Link:
        % http://books.google.com/books?id=aqjU_NHyre4C&lpg=PP1&dq=Global%20
        % earth%20physics%3A%20a%20handbook%20of%20physical%20constants&pg=PA49
        % #v=onepage&q&f=false
        % (except equation 6 is missing a required 1/sin(theta) and m; correct
        % equations on page 5 (equations 3a-3c) of:
        % http://hanspeterschaub.info/Papers/UnderGradStudents/
        % MagneticField.pdf)
        if m == 0 % Implies h = 0, so only coefficient in gh is g
            coef = a_r*gh(coefindex); %*cos(0*phi) = 1
            Br = Br + (n+1)*coef*P(Pindex);
            Bt = Bt - coef*dP(Pindex);
            % Bp is 0 for m = 0.
            coefindex = coefindex + 1; % Only need to skip over g this time.
        else
            coef = a_r*(gh(coefindex)*cosphi(m) + gh(coefindex+1)*sinphi(m));
            Br = Br + (n+1)*coef*P(Pindex);
            Bt = Bt - coef*dP(Pindex);
            if sintheta == 0 % Use different formula when dividing by 0.
                Bp = Bp - costheta*a_r*(-gh(coefindex)*sinphi(m) + ...
                    gh(coefindex+1)*cosphi(m))*dP(Pindex);
            else
                Bp = Bp - 1/sintheta*a_r*m*(-gh(coefindex)*sinphi(m) + ...
                    gh(coefindex+1)*cosphi(m))*P(Pindex);
            end
            coefindex = coefindex + 2; % Skip over g and h this time.
        end

        % Increment m.
        m = m + 1;

    end

    % Convert from spherical to (x,y,z) = (North,East,Down).
    Bx = -Bt;
    By = Bp;
    Bz = -Br;

    % Convert back to geodetic coordinates
    Bx_old = Bx;
    Bx = Bx.*cd + Bz.*sd;
    Bz = Bz.*cd - Bx_old.*sd;
    
    % Build output vector NED
    B = [Bx;By;Bz];
    
    % convert to ecef from ned
    Ned2ecef = frameRotationNed2ecef(sinlat, coslat, sinphi(1), cosphi(1));
    Becef = Ned2ecef*B;
end
