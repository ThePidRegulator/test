% constants.m
%
% Created: 28.01.2015 16:56:50
% Author: Antoine Pignede
%
% This script shall define ALL the global constants needed anywhere in the
%   prediction algorithms for sun, orbital position and geomagnetic field.

    load tle;

    global DEG2RAD
    DEG2RAD = pi / 180.0;

%% Excerpt of constmath.m from Vallado (2013) Fundamentals of Astrodynamics and Applications
    global RAD2DEG TWOPI
    RAD2DEG = 180.0 / pi;
    TWOPI = 2.0 * pi;
    %HALFPI = pi * 0.5;

%% Excerpt of constastro.m from Vallado (2013) Fundamentals of Astrodynamics and Applications
    % WGS-84/EGM-96 constants used here
    global REQU FLAT RPOL RMEAN MU
    REQU       = 6378.137;         % km semimajor axis
    FLAT       = 1.0/298.257223563;
    RPOL       = REQU*(1 - FLAT);  % km semiminor axis
    RMEAN      = 6371.2;           % km mean earth radius
    %OMEGAEARTH = 7.292115e-5;     % rad/s
    MU         = 398600.4418;      % km3/s2

    % derived constants from the base values
    global ECCEARTH ECCEARTHSQRD TUSEC TUMIN VELKMPS
    ECCEARTH = sqrt(2.0*FLAT - FLAT^2);
    ECCEARTHSQRD = ECCEARTH^2;
    TUSEC = sqrt(REQU^3/MU);
    TUMIN = TUSEC / 60.0;
    %TUDAY = TUSEC / 86400.0;
    %OMEGAEARTHRADPTU = OMEGAEARTH * TUSEC;
    %OMEGAEARTRADPMIN = OMEGAEARTH * 60.0;
    VELKMPS = sqrt(MU / REQU);
    %VELRADPMIN = VELKMPS * 60.0/RE;
    %for afspc
    %velkmps1 = velradpmin*6378.135/60.0   7.90537051051763
    %mu1 = velkmps*velkmps*6378.135        3.986003602567418e+005        
    %DEGPSEC = (180.0 / pi) / TUSEC;
    %RADPDAY = 2.0 * pi * 1.002737909350795;

%% Excerpt fo getgravc.m from Vallado (2013) Fundamentals of Astrodynamics and Applications
    % ------------ wgs-84 constants ------------
    global XKE J2 J3 J4 J3OJ2
    XKE    = 1.0 / TUMIN;
    J2     =   0.00108262998905;
    J3     =  -0.00000253215306;
    J4     =  -0.00000161098761;
    J3OJ2  =  J3 / J2;

%% Excerpt of iau80in.m from Vallado (2013) Fundamentals of Astrodynamics and Applications
    global IAR80 RAR80
    % 0.0001" to rad
    convrt = 0.0001 * pi / (180*3600.0);

    load nut80.dat;
    IAR80 = nut80(:,1:5);
    RAR80 = nut80(:,6:9);

    for i=1:106
       for j=1:4
           RAR80(i,j)= RAR80(i,j) * convrt;
       end
    end

%% Polar motion estimation coefficients
    % Values taken from IERS Bulletin A 12 February 2015 Vol. XXVIII No. 007 http://www.iers.org/IERS/EN/Publications/Bulletins/bulletins.html
    global MJDPOLAR XPOLAR YPOLAR
    MJDPOLAR = 57065;
    XPOLAR = [0.1095;-0.0737;-0.0326;-0.0373;0.0080];
    YPOLAR = [0.3519;-0.0293;0.0661;0.0080;0.0373];

%% IGRF12 coefficients. Inspired by loadigrfcoefs.m from Compston (2011) Matlab File Exchange International Geomagnetic Reference Field (IGRF) Model
    global IGRFGH IGRFGHSV
    igrfdata = xlsread('IGRF12coeffsMatlabFormat.xls','igrf12coeffsMatlabFormat');
    IGRFGH = igrfdata(:,3);
    IGRFGHSV = igrfdata(:,4);

%% WMM2015 coefficients
    global WMMG WMMGSV WMMH WMMHSV
    wmmdata = dlmread('WMM2015coeffsMatlabFormat.COF');
    WMMG = wmmdata(:,3);
    WMMH = wmmdata(:,4);
    WMMGSV = wmmdata(:,5);
    WMMHSV = wmmdata(:,6);
