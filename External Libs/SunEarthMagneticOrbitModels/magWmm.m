%Geomagnetic Vector Prediction
% function magWmm
%
% Created: 07.02.2015 08:46:08
% Author: Antoine Pignede
%
% Calculates the components of the Earth's magnetic field using the
%   World Magnetic Model (WMM).
%   Note that the other parameters the WMM gives can be computed from BX,
%   BY, and BZ as:
%     Horizonal intensity: hypot(BX, BY) (i.e., sqrt(BX.^2 + BY.^2) )
%     Total intensity: hypot(BX, hypot(BY, BZ))
%     Declination: atan2(BY, BX)
%     Inclination: atan(BZ./hypot(BX, BY))
%
% Copied and edited from the official C source code available at http://www.ngdc.noaa.gov/geomag/WMM/soft.shtml
% 
% Inputs:
%    reci        - position vector eci            km
%    timeUTC     - struct containing:
%      year        - year                           1900 .. 2100
%      mon         - month                          1 .. 12
%      day         - day                            1 .. 28,29,30,31
%      hr          - universal time hour            0 .. 23
%      min         - universal time min             0 .. 59
%      sec         - universal time sec             0.0 .. 59.999
% 
% Outputs:
%    Beci        - geomagnetic field vector in nanotesla (nT) in eci frame
%
% Coupling
%    timeDatetime2years     - convert date and time values to fractional years
%    frameEci2ecef     - transform eci position to ecef position
%    frameEcef2geod    - transform ecef position to geodetic position
%
% Note that the original functions have been merged into this wmm function.
% 
% See also
%    magIgrf.m  - alternative geomagnetic vector predictor

function [Beci] = magWmm(reci,timeUTC)
global DEG2RAD RAD2DEG REQU RMEAN ECCEARTHSQRD WMMG WMMGSV WMMH WMMHSV

    % Convert time to fractional years.
    time = timeDatetime2years(timeUTC);

    % Transform eci inputvector to Geodetic coordinates
    recef = frameEci2ecef(reci,timeUTC);
    rgeod = frameEcef2geod(recef);
    latitude = rgeod(1);
    longitude = rgeod(2);
    altitude = rgeod(3);
    
    % Convert Geodetic coordinates to Spherical coordinates    
    coslat = cos(DEG2RAD*latitude);
    sinlat = sin(DEG2RAD*latitude);
    % compute the local radius of curvature on the WGS-84 reference ellipsoid
    rc = REQU / sqrt(1.0 - ECCEARTHSQRD * sinlat * sinlat);
    % compute ECEF Cartesian coordinates of specified point (for longitude=0)
    xp = (rc + altitude) * coslat;
    zp = (rc * (1.0 - ECCEARTHSQRD) + altitude) * sinlat;
    % compute spherical radius and angle lambda and phi of specified point
    r = sqrt(xp * xp + zp * zp);
    phi = RAD2DEG*asin(zp / r);
    lambda = longitude;
    coslambda = cos(DEG2RAD * lambda);
    sinlambda = sin(DEG2RAD * lambda);
    cosphi = sin(DEG2RAD * phi);
    sinphi = sin(DEG2RAD * phi);
    
    % Time change the Model coefficients from the base year of the model using secular variation coefficients
    g = WMMG + WMMGSV*(time - 2015);
    h = WMMH + WMMHSV*(time - 2015);
    nMax = 12;
    NumTerms = (nMax + 1) * (nMax + 2) / 2;
    
    % Compute Spherical variables: (a/r)^(n+2), cos_m(lamda) and sin_m(lambda) for spherical harmonic summations
    RelativeRadiusPower = zeros(nMax+1,1);
    RelativeRadiusPower(1) = (RMEAN / r) * (RMEAN / r);
    for n = 1:nMax
        RelativeRadiusPower(n+1) = RelativeRadiusPower(n) * (RMEAN / r);
    end
    cosmlambda = zeros(nMax+1,1);
    sinmlambda = zeros(nMax+1,1);
    cosmlambda(1) = 1.0;
    sinmlambda(1) = 0.0;
    cosmlambda(2) = coslambda;
    sinmlambda(2) = sinlambda;
    for m = 2:nMax
        cosmlambda(m+1) = cosmlambda(m) * coslambda - sinmlambda(m) * sinlambda;
        sinmlambda(m+1) = cosmlambda(m) * sinlambda + sinmlambda(m) * coslambda;
    end
    
    % Compute all of the Schmidt-semi normalized associated Legendre functions up to degree nMax
    Pcup = zeros(NumTerms,1);
    dPcup = zeros(NumTerms,1);
    Pcup(1) = 1.0;
    dPcup(1) = 0.0;
    z = sqrt((1.0 - sinphi) * (1.0 + sinphi));
    % First, Compute the Gauss-normalized associated Legendre functions
    for n = 1:nMax
        for m = 0:n
            index = n * (n + 1) / 2 + m;
            if(n == m)
                index1 = (n - 1) * n / 2 + m - 1;
                Pcup(index+1) = z * Pcup(index1+1);
                dPcup(index+1) = z * dPcup(index1+1) + sinphi * Pcup(index1+1);
            elseif((n == 1) && (m == 0))
                index1 = (n - 1) * n / 2 + m;
                Pcup(index+1) = sinphi * Pcup(index1+1);
                dPcup(index+1) = sinphi * dPcup(index1+1) - z * Pcup(index1+1);
            elseif((n > 1) && (n ~= m))
                index1 = (n - 2) * (n - 1) / 2 + m;
                index2 = (n - 1) * n / 2 + m;
                if(m > (n - 2))
                    Pcup(index+1) = sinphi * Pcup(index2+1);
                    dPcup(index+1) = sinphi * dPcup(index2+1) - z * Pcup(index2+1);
                else
                    k = (((n - 1) * (n - 1)) - (m * m)) / ((2 * n - 1) * (2 * n - 3));
                    Pcup(index+1) = sinphi * Pcup(index2+1) - k * Pcup(index1+1);
                    dPcup(index+1) = sinphi * dPcup(index2+1) - z * Pcup(index2+1) - k * dPcup(index1+1);
                end
            end
        end
    end
    % Compute the ration between the Schmidt quasi-normalized associated Legendre functions and the Gauss-normalized version
    schmidtQuasiNorm = zeros(NumTerms,1);
    schmidtQuasiNorm(1) = 1.0;
    for n = 1:nMax
        index = n * (n + 1) / 2;
        index1 = (n - 1) * n / 2;
        schmidtQuasiNorm(index+1) = schmidtQuasiNorm(index1+1) * (2 * n - 1) / n;
        for m = 1:n
            index = n * (n + 1) / 2 + m;
            index1 = n * (n + 1) / 2 + m - 1;
            oneOrTwo = 1;
            if m == 1
                oneOrTwo = 2;
            end
            schmidtQuasiNorm(index+1) = schmidtQuasiNorm(index1+1) * sqrt(((n - m + 1) * (oneOrTwo)) /  (n + m));
        end
    end
    % Convert the Gauss-normalized associated Legendre functions to the Schmidt quasi-normalized version
    for n = 1:nMax
        for m = 0:n
            index = n * (n + 1) / 2 + m;
            Pcup(index+1) = Pcup(index+1) * schmidtQuasiNorm(index+1);
            dPcup(index+1) = -dPcup(index+1) * schmidtQuasiNorm(index+1);
        end
    end
    
    % Compute Geomagnetic Field Elements X, Y and Z in Spherical coordinate system using spherical harmonic summation
    Bx = 0;
    By = 0;
    Bz = 0;
    for n = 1:nMax
        for m = 0:n
            index = n * (n + 1) / 2 + m;
            Bz = Bz - (RelativeRadiusPower(n+1) * (g(index) * cosmlambda(m+1) + h(index) * sinmlambda(m+1)) * (n + 1) * Pcup(index+1));
            By = By + (RelativeRadiusPower(n+1) * (g(index) * sinmlambda(m+1) - h(index) * cosmlambda(m+1)) * m * Pcup(index+1));
            Bx = Bx - (RelativeRadiusPower(n+1) * (g(index) * cosmlambda(m+1) + h(index) * sinmlambda(m+1)) * dPcup(index+1));
        end
    end
    if(abs(cosphi) > 1.0e-10)
        By = By / cosphi;
    else % Special calculation for component By at Geographic poles        
        PcupS = zeros(nMax+1,1);
        PcupS(1) = 1;
        schmidtQuasiNorm1 = 1.0;
        By = 0.0;
        % Compute the ratio between the Gauss-normalized associated Legendre functions and the Schmidt quasi-normalized version
        for n = 1:nMax
            index = n * (n + 1) / 2 + 1;
            schmidtQuasiNorm2 = schmidtQuasiNorm1 * (2 * n - 1) / n;
            schmidtQuasiNorm3 = schmidtQuasiNorm2 * sqrt((n * 2) / (n + 1));
            schmidtQuasiNorm1 = schmidtQuasiNorm2;
            if(n == 1)
                PcupS(n+1) = PcupS(n);
            else
                k = (((n - 1) * (n - 1)) - 1) / ((2 * n - 1) * (2 * n - 3));
                PcupS(n+1) = sinphi * PcupS(n) - k * PcupS(n-1);
            end
            By = By + (RelativeRadiusPower(n+1) * (g(index) * sinmlambda(2) - h(index) * cosmlambda(2)) * PcupS(n+1) * schmidtQuasiNorm3);
        end       
    end
    
    % Convert back to geodetic coordinates
    sin_psi = sin(DEG2RAD * (phi - latitude));
    cos_psi = cos(DEG2RAD * (phi - latitude));
    Bx_old = Bx;
    Bx = Bx_old * cos_psi - Bz * sin_psi;
    Bz = Bx_old * sin_psi + Bz * cos_psi;

    % Build output vector NED
    B = [Bx;By;Bz];
    
    % Convert to ECI from Fossen (2011) Handbook of Marine Craft Hydrodynamics and Motion Control
    Ned2ecef = frameRotationNed2ecef(sinlat, coslat, sinlambda(1), coslambda(1));
    Becef = Ned2ecef*B;
    Beci = frameEcef2eci(Becef,timeUTC);
end
