% testMatlab.m
%
% Created: 10.02.2015 13:10:11
% Author: Antoine Pignede
%
% Use this script to test the Matlab functions. Each functions is used at 
%   least once with known input and output values given in the sources.

    clear all;
    close all;
    testCount = 0;
    errorCount = 0;
    absTol = 0.001; % tolerance for the absolute difference when comparing matrices
    relTol = 0.01; % tolerance for the relative difference when comparing matrices

%% Define all global constants
    constants;
    
%% Test the Time Conversions
    % Excerpt of ex3_1415.m from Vallado (2013) Fundamentals of Astrodynamics and Applications
    timeSol.year=2004;
    timeSol.mon = 4;
    timeSol.day = 6;
    timeSol.hr =  7;
    timeSol.min= 51;
    timeSol.sec= 28.386009;
    jd = timeDatetime2jd(timeSol); % TEST!
    testCount = testCount + 1;
    jdSol = 2453101.8274118751;
    if abs(jd - jdSol) > (1/86400) % If difference greater than 1s
        errorCount = errorCount +1;
        fprintf(1,'timeDatetime2jd error! %18.10f not %18.10f\n',jd,jdSol);
    else
        fprintf(1,'Test %d timeDatetime2jd OK!\n',testCount);
    end
    
    time = timeJd2datetime(jdSol); % TEST!
    testCount = testCount + 1;
    if (time.year ~= timeSol.year) || (time.mon ~= timeSol.mon) || (time.day ~= timeSol.day) || (time.hr ~= timeSol.hr) || (time.min ~= timeSol.min) || (abs(time.sec - timeSol.sec) > 0.1)
        errorCount = errorCount +1;
        fprintf(1,'timeJd2datetime error! time not timeSol\n');
        time
        timeSol
    else
        fprintf(1,'Test %d timeJd2datetime OK!\n',testCount);
    end
    
    days = timeDatetime2days(timeSol); % TEST!
    testCount = testCount + 1;
    daysSol = 97.3274118751;
    if abs(days - daysSol) > (1/86400) % If difference greater than 1s
        errorCount = errorCount +1;
        fprintf(1,'timeDatetime2days error! %18.10f not %18.10f\n',days,daysSol);
    else
        fprintf(1,'Test %d timeDatetime2days OK!\n',testCount);
    end
    
    time = timeDays2datetime(timeSol.year,daysSol); % TEST
    testCount = testCount + 1;
    if (time.year ~= timeSol.year) || (time.mon ~= timeSol.mon) || (time.day ~= timeSol.day) || (time.hr ~= timeSol.hr) || (time.min ~= timeSol.min) || (abs(time.sec - timeSol.sec) > 1)
        errorCount = errorCount +1;
        fprintf(1,'timeDays2datetime error! time not timeSol\n');
        time
        timeSol
    else
        fprintf(1,'Test %d timeDays2datetime OK!\n',testCount);
    end
    
    years = timeDatetime2years(timeSol); % TEST!
	testCount = testCount + 1;
	yearsSol =  2004.2631896499;
	if abs(years - yearsSol) > (1.0/86400) % If difference greater than 1s
		errorCount = errorCount +1;
		fprintf(1,'timeDatetime2years error! %18.10f not %18.10f\n',years,yearsSol);
    else
        fprintf(1,'Test %d timeDatetime2years OK!\n',testCount);
    end

	time = timeYears2datetime(yearsSol); % TEST
	testCount = testCount + 1;
	if (time.year ~= timeSol.year) || (time.mon ~= timeSol.mon) || (time.day ~= timeSol.day) || (time.hr ~= timeSol.hr) || (time.min ~= timeSol.min) || (abs(time.sec - timeSol.sec) > 1)
        errorCount = errorCount +1;
        fprintf(1,'timeYears2datetime error! time not timeSol\n');
        time
        timeSol
    else
        fprintf(1,'Test %d timeYears2datetime OK!\n',testCount);
    end
    
    jc = timeJd2jc(jdSol); % TEST
    testCount = testCount + 1;
    jcSol = 0.0426236116;
    if abs(jc - jcSol) > (1/3.1558e+09) % If difference greater than 1s
        errorCount = errorCount +1;
        fprintf(1,'timeJd2jc error! %18.10f not %18.10f\n',jc,jcSol);
    else
        fprintf(1,'Test %d timeJd2jc OK!\n',testCount);
    end
    
    jd = timeJc2jd(jc); % TEST
    testCount = testCount + 1;
    if abs(jd - jdSol) > (1/86400) % If difference greater than 1s
        errorCount = errorCount +1;
        fprintf(1,'timeJc2jd error! %18.10f not %18.10f\n',jd,jdSol);
    else
        fprintf(1,'Test %d timeJc2jd OK!\n',testCount);
    end
    
    timeTT = timeUtc2tt(timeSol); % TEST
    testCount = testCount + 1;
    timeTTSol.year=2004;
    timeTTSol.mon = 4;
    timeTTSol.day = 6;
    timeTTSol.hr =  7;
    timeTTSol.min= 52;
    timeTTSol.sec= 32.570009;
    if (timeTT.year ~= timeTTSol.year) || (timeTT.mon ~= timeTTSol.mon) || (timeTT.day ~= timeTTSol.day) || (timeTT.hr ~= timeTTSol.hr) || (timeTT.min ~= timeTTSol.min) || (abs(timeTT.sec - timeTTSol.sec) > 1)
        errorCount = errorCount +1;
        fprintf(1,'timeUtc2tt error! timeTT not timeTTSol\n');
        timeTT
        timeTTSol
    else
        fprintf(1,'Test %d timeUtc2tt OK!\n',testCount);
    end
    
    time = timeTt2utc(timeTTSol); % TEST
    testCount = testCount + 1;
    if (time.year ~= timeSol.year) || (time.mon ~= timeSol.mon) || (time.day ~= timeSol.day) || (time.hr ~= timeSol.hr) || (time.min ~= timeSol.min) || (abs(time.sec - timeSol.sec) > 1)
        errorCount = errorCount +1;
        fprintf(1,'timeTt2utc error! time not timeSol\n');
        time
        timeSol
    else
        fprintf(1,'Test %d timeTt2utc OK!\n',testCount);
    end
    
    gst = timeGstime(jdSol); % TEST
	testCount = testCount + 1;
	gstSol = 5.4594896669;
	if (abs(gst - gstSol) > absTol)
		errorCount = errorCount +1;
		fprintf(1,'timeGstime error! time not timeSol\n');
        gst
        gstSol
    else
        fprintf(1,'Test %d timeGstime OK!\n\n',testCount);
    end
    
%% Test the Frame Transformations
    tttSol = 0.042623631889;
	[l,l1,f,d,omega] = frameFundarg(tttSol); % TEST
	testCount = testCount + 1;
	lSol = 5.4962562692;
	l1Sol = 1.6046448385;
	fSol = 2.9512994807;
	dSol = 3.4339639764;
	omegaSol = 0.7435913620;
	if ((abs(l - lSol) > absTol) || (abs(l1 - l1Sol) > absTol) || (abs(f - fSol) > absTol) || (abs(d - dSol) > absTol) || (abs(omega - omegaSol) > absTol))
		errorCount = errorCount +1;
		fprintf(1,'frameFundarg error! delaunay not delaunaySol\n');
		l
        lSol
        l1
        l1Sol
        f
        fSol
        d
        dSol
        omega
        omegaSol
    else
        fprintf('Test %d frameFundarg OK!\n',testCount);
    end
    
    % Excerpt of ex3_1415.m from Vallado (2013) Fundamentals of Astrodynamics and Applications
    prec = framePrecess(tttSol); % TEST
    testCount = testCount + 1;
    precSol = [0.99999945998 0.00095314992 0.00041417739;-0.00095314992 0.99999954575 -0.00000019739;-0.00041417739 -0.00000019739 0.99999991423];
    absError = prec(:)-precSol(:);
    relError = absError./prec(:);
    relError(~isfinite(relError)) = 0;   % Sets Inf and NaN to 0
    same = any( (abs(absError) > absTol) | (abs(relError) > relTol) );
    if same
        errorCount = errorCount +1;
        fprintf(1,'framePrecess error! prec not precSol\n');
        prec
        precSol
    else
        fprintf(1,'Test %d framePrecess OK!\n',testCount);
    end
    
    [deltapsi,meaneps,omega,nut] = frameNutation(tttSol); % TEST
    testCount = testCount + 1;
    nutSol = [0.99999999823 -0.00005461818 -0.00002367925;0.00005461734 0.99999999788 0.00003545941;0.00002368119 -0.00003545812 0.99999999909];
    omegaSol = 0.7435907904;
    meanepsSol = 0.4090831301;
    deltapsiSol = -0.0000597833;
    absError = nut(:)-nutSol(:);
    relError = absError./nut(:);
    relError(~isfinite(relError)) = 0;   % Sets Inf and NaN to 0
    same = any( (abs(absError) > absTol) | (abs(relError) > relTol) );
    if same
        errorCount = errorCount +1;
        fprintf(1,'frameNutation error! nut not nutSol\n');
        nut
        nutSol
        fprintf(1,'omega=%18.10f omegaSol=%18.10f\n',omega,omegaSol);
        fprintf(1,'meaneps=%18.10f meanepsSol=%18.10f\n',meaneps,meanepsSol);
        fprintf(1,'deltapsi=%18.10f deltapsiSol=%18.10f\n',deltapsi,deltapsiSol);
    else
        fprintf(1,'Test %d frameNutation OK!\n',testCount);
    end
    
    if  (abs(omega - omegaSol) > 0.001) || (abs(meaneps - meanepsSol) > 0.001) || (abs(deltapsi - deltapsiSol) > 0.001)
        fprintf(1,'frameNutation warning!\n');
        nut
        nutSol
        fprintf(1,'omega=%18.10f omegaSol=%18.10f\n',omega,omegaSol);
        fprintf(1,'meaneps=%18.10f meanepsSol=%18.10f\n',meaneps,meanepsSol);
        fprintf(1,'deltapsi=%18.10f deltapsiSol=%18.10f\n',deltapsi,deltapsiSol);
    end
    
    st  = frameSidereal(jdSol,deltapsiSol,meanepsSol,omegaSol); % TEST
    testCount = testCount + 1;
    stSol = [0.67952793705 0.73364963216 0.00000000000;-0.73364963216 0.67952793705 0.00000000000;0.00000000000 0.00000000000 1.00000000000];
    absError = st(:)-stSol(:);
    relError = absError./st(:);
    relError(~isfinite(relError)) = 0;   % Sets Inf and NaN to 0
    same = any( (abs(absError) > absTol) | (abs(relError) > relTol) );
    if same
        errorCount = errorCount +1;
        fprintf(1,'frameSidereal error! st not stSol\n');
        st
        stSol
    else
        fprintf(1,'Test %d frameSidereal OK!\n',testCount);
    end
    
    pm = framePolarm(jdSol); % TEST
    testCount = testCount + 1;
    pmSol = [1.00000000000 0.00000000000 0.00000068205;-0.00000000000 1.00000000000 0.00000161593;-0.00000068205 -0.00000161593 1.00000000000];
    absError = pm(:)-pmSol(:);
    relError = absError./pm(:);
    relError(~isfinite(relError)) = 0;   % Sets Inf and NaN to 0
    same = any( (abs(absError) > absTol) | (abs(relError) > relTol) );
    if same
        errorCount = errorCount +1;
        fprintf(1,'framePolarm error! pm not pmSol\n');
        pm
        pmSol
    else
        fprintf(1,'Test %d framePolarm OK!\n',testCount);
    end
    
    nutteme = frameTruemean(tttSol); % TEST
    testCount = testCount + 1;
    nuttemeSol = [1.00000000000 -0.00000000897 -0.00000041328;0.00000000897 1.00000000000 0.00000061887;0.00000041328 -0.00000061887 1.00000000000];
    absError = nutteme(:)-nuttemeSol(:);
    relError = absError./nutteme(:);
    relError(~isfinite(relError)) = 0;   % Sets Inf and NaN to 0
    same = any( (abs(absError) > absTol) | (abs(relError) > relTol) );
    if same
        errorCount = errorCount +1;
        fprintf(1,'frameTruemean error! nutteme not nuttemeSol\n');
        nutteme
        nuttemeSol
    else
        fprintf(1,'Test %d frameTruemean OK!\n',testCount);
    end
    
    recefSol = [-1033.4793830;7901.2952754;6380.3565958];
    reci = frameEcef2eci(recefSol,timeSol); % TEST
    testCount = testCount + 1;
    reciSol = [5102.5089579;6123.0114007;6378.1369282];
    absError = reci(:)-reciSol(:);
    relError = absError./reci(:);
    relError(~isfinite(relError)) = 0;   % Sets Inf and NaN to 0
    same = any( (abs(absError) > absTol) | (abs(relError) > relTol) );
    if same
        errorCount = errorCount +1;
        fprintf(1,'frameEcef2eci error! reci not reciSol\n');
        reci
        reciSol
    else
        fprintf(1,'Test %d frameEcef2eci OK!\n',testCount);
    end
    
    recef = frameEci2ecef(reciSol,timeSol); % TEST
    testCount = testCount + 1;
    absError = recef(:)-recefSol(:);
    relError = absError./recef(:);
    relError(~isfinite(relError)) = 0;   % Sets Inf and NaN to 0
    same = any( (abs(absError) > absTol) | (abs(relError) > relTol) );
    if same
        errorCount = errorCount +1;
        fprintf(1,'frameEci2ecef error! recef not recefSol\n');
        recef
        recefSol
    else
        fprintf(1,'Test %d frameEci2ecef OK!\n',testCount);
    end
    
    rteme = frameEci2teme(reciSol,timeSol); % TEST
    testCount = testCount + 1;
    rtemeSol = [5094.1801621;6127.6446595;6380.3445327];
    absError = rteme(:)-rtemeSol(:);
    relError = absError./rteme(:);
    relError(~isfinite(relError)) = 0;   % Sets Inf and NaN to 0
    same = any( (abs(absError) > absTol) | (abs(relError) > relTol) );
    if same
        errorCount = errorCount +1;
        fprintf(1,'frameEci2teme error! rteme not rtemeSol\n');
        rteme
        rtemeSol
    else
        fprintf(1,'Test %d frameEci2teme OK!\n',testCount);
    end
    
    reci = frameTeme2eci(rtemeSol,timeSol); % TEST
    testCount = testCount + 1;
    absError = reci(:)-reciSol(:);
    relError = absError./reci(:);
    relError(~isfinite(relError)) = 0;   % Sets Inf and NaN to 0
    same = any( (abs(absError) > absTol) | (abs(relError) > relTol) );
    if same
        errorCount = errorCount +1;
        fprintf(1,'frameTeme2eci error! reci not reciSol\n');
        reci
        reciSol
    else
        fprintf(1,'Test %d frameTeme2eci OK!\n',testCount);
    end
    
   	% Excerpt of ex3_3.m from Vallado (2013) Fundamentals of Astrodynamics and Applications
	recefSol = [6524.834;6862.875;6448.296];
	rgeod = frameEcef2geod(recefSol); % TEST
	testCount = testCount + 1;
	rgeodSol = [34.3524952;46.4464169;5085.2187362];
    absError = rgeod(:)-rgeodSol(:);
    relError = absError./rgeod(:);
    relError(~isfinite(relError)) = 0;   % Sets Inf and NaN to 0
	same = any( (abs(absError) > absTol) | (abs(relError) > relTol) );
    if same
        errorCount = errorCount +1;
        fprintf(1,'frameEcef2geod error! rgeod not rgeodSol\n');
        rgeod
        rgeodSol
    else
        fprintf(1,'Test %d frameEcef2geod OK!\n',testCount);
    end

	recef = frameGeod2ecef(rgeodSol); % TEST
	testCount = testCount + 1;
	absError = recef(:)-recefSol(:);
    relError = absError./recef(:);
    relError(~isfinite(relError)) = 0;   % Sets Inf and NaN to 0
	same = any( (abs(absError) > absTol) | (abs(relError) > relTol) );
    if same
        errorCount = errorCount +1;
        fprintf(1,'frameGeod2ecef error! recef not recefSol\n');
        recef
        recefSol
    else
        fprintf(1,'Test %d frameGeod2ecef OK!\n',testCount);
    end
    
    % Verify frameRotationNed2ecef with the Marine Systems Simulator (MSS) toolbox
    sinlat = sin(rgeodSol(1));
    coslat = cos(rgeodSol(1));
    sinlon = sin(rgeodSol(2));
    coslon = cos(rgeodSol(2));
    Ned2ecef =  frameRotationNed2ecef(sinlat, coslat, sinlon, coslon); % TEST
    testCount = testCount + 1;
    Ned2ecefSol = [0.1586295368 -0.6268260173 -0.7628406217;-0.1276159144 -0.7791592546 0.6136978363;-0.9790561008 0.0000000000 -0.2035906470];
    absError = Ned2ecef(:)-Ned2ecefSol(:);
    relError = absError./Ned2ecef(:);
    relError(~isfinite(relError)) = 0;   % Sets Inf and NaN to 0
    same = any( (abs(absError) > absTol) | (abs(relError) > relTol) );
    if same
        errorCount = errorCount +1;
        fprintf(1,'frameRotationNed2ecef error! Ned2ecef not Ned2ecefSol\n');
        Ned2ecef
        Ned2ecefSol
    else
        fprintf(1,'Test %d frameRotationNed2ecef OK!\n\n',testCount);
    end    
    
%% Test the Sun Vector Prediction
    % Excerpt of ex5_1.m from Vallado (2013) Fundamentals of Astrodynamics and Applications
    timeSol.year=2006;
	timeSol.mon = 4;
	timeSol.day = 2;
	timeSol.hr =  0;
	timeSol.min= 0;
	timeSol.sec= 0;
    rsun = sun(timeSol); % TEST
    testCount = testCount + 1;
    rsunSol = [0.9780491365;0.1911814374;0.0828827162];
    rsunSol = rsunSol/norm(rsunSol);
    absError = rsun(:)-rsunSol(:);
    relError = absError./rsun(:);
    relError(~isfinite(relError)) = 0;   % Sets Inf and NaN to 0
    same = any( (abs(absError) > absTol) | (abs(relError) > relTol) );
    if same
        errorCount = errorCount +1;
        fprintf(1,'sun error! rsun not rsunSol\n');
        rsun
        rsunSol
    else
        fprintf(1,'Test %d sun OK!\n\n',testCount);
    end
    
%% Test the Orbit Position Prediction
    % Original STR#3 SGP4 test from ValladoEtAl (2006) Revisiting Spacetrack Report no3
    tleSol = '1 88888U          80275.98708465  .00073094  13844-3  66816-4 0    872 88888  72.8435 115.9689 0086731  52.6988 110.5714 16.05824518  1058';
    satrec = orbitTwoline2rv(tleSol);
    fprintf(1,'satrec.error after initialization: "%s"\n',satrec.error);
    timeSol.year=1980;
    timeSol.mon = 10;
    timeSol.day = 2;
    timeSol.hr =  1;
    timeSol.min= 41;
    timeSol.sec= 24.113771;
    [satrec, reci, veci] = orbitSgp4(satrec,timeSol); % TEST  % TEST
    fprintf(1,'satrec.error after SGP4: "%s"\n',satrec.error);
    testCount = testCount + 1;
    reciSol = [1022.4233730086;2290.9626955932;-6189.6432266006];
    absError = reci(:)-reciSol(:);
    relError = absError./reci(:);
    relError(~isfinite(relError)) = 0;   % Sets Inf and NaN to 0
    same = any( (abs(absError) > absTol) | (abs(relError) > relTol) );
    if same
        errorCount = errorCount +1;
        fprintf(1,'orbitSgp4 position error! reci not reciSol\n');
        reci
        reciSol
    else
        fprintf(1,'Test %d orbitSgp4 position OK!\n',testCount);
    end
    testCount = testCount +1;
    veciSol = [-3.7777586815;6.4513394443;1.8209505936];
    absError = veci(:)-veciSol(:);
    relError = absError./veci(:);
    relError(~isfinite(relError)) = 0;   % Sets Inf and NaN to 0
    same = any( (abs(absError) > absTol) | (abs(relError) > relTol) );
    if same
        errorCount = errorCount +1;
        fprintf(1,'orbitSgp4 velocity error! veci not veciSol\n');
        veci
        veciSol
    else
        fprintf(1,'Test %d orbitSgp4 velocity OK!\n',testCount);
    end
    
    % TEME example from ValladoEtAl (2006) Revisiting Spacetrack Report no3
	tleSol2 = '1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  47532 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667';
	satrec = orbitTwoline2rv(tleSol2);
	fprintf(1,'satrec.error after initialization: "%s"\n',satrec.error);
	timeSol.year=2000;
	timeSol.mon = 6;
	timeSol.day = 30;
	timeSol.hr =  18;
	timeSol.min= 50;
	timeSol.sec= 19.733571;
	[satrec, reci, veci] = orbitSgp4(satrec,timeSol); % TEST  % TEST
	fprintf(1,'satrec.error after SGP4: "%s"\n',satrec.error);
	testCount = testCount + 1;
	reciSol = [-9059.9413786;4659.6972000;813.9588875];
	absError = reci(:)-reciSol(:);
    relError = absError./reci(:);
    relError(~isfinite(relError)) = 0;   % Sets Inf and NaN to 0
    same = any( (abs(absError) > absTol) | (abs(relError) > relTol) );
    if same
        errorCount = errorCount +1;
        fprintf(1,'orbitSgp4 position error! reci not reciSol\n');
        reci
        reciSol
    else
        fprintf(1,'Test %d orbitSgp4 position OK!\n',testCount);
    end
    testCount = testCount +1;
    veciSol = [-2.233348094;-4.110136162;-3.157394074];
    absError = veci(:)-veciSol(:);
    relError = absError./veci(:);
    relError(~isfinite(relError)) = 0;   % Sets Inf and NaN to 0
    same = any( (abs(absError) > absTol) | (abs(relError) > relTol) );
    if same
        errorCount = errorCount +1;
        fprintf(1,'orbitSgp4 velocity error! veci not veciSol\n');
        veci
        veciSol
    else
        fprintf(1,'Test %d orbitSgp4 velocity OK!\n\n',testCount);
    end
    
%% Test the Geomagnetic Vector Prediction
    % Official WMM 2015 test values
    timeSol.year=2015;
    timeSol.mon = 1;
    timeSol.day = 1;
    timeSol.hr =  0;
    timeSol.min= 0;
    timeSol.sec= 0.0;
    reciSol = [-186.3979236310;1093.5121033253;6259.8764001403];
    tic;
    Beci = magWmm(reciSol,timeSol); % TEST
    toc;
    testCount = testCount + 1;
    BeciSol = [3174.4799983588;-15647.9734763282;-52460.0435315694];
    absError = Beci(:)-BeciSol(:);
    relError = absError./Beci(:);
    relError(~isfinite(relError)) = 0;   % Sets Inf and NaN to 0
    same = any( (abs(absError) > absTol) | (abs(relError) > relTol) );
    if same
        errorCount = errorCount +1;
        fprintf(1,'magWmm error! Beci not BeciSol\n');
        Beci
        BeciSol
    else
        fprintf(1,'Test %d magWmm OK!\n',testCount);
    end
    
    % Magnetic Field Calculator http://www.ngdc.noaa.gov/geomag-web/?model=igrf#igrfwmm
    tic;
    Beci = magIgrf(reciSol,timeSol); % TEST
    toc;
    testCount = testCount + 1;
    BeciSol = [3176.4834170946;-15651.8044943423;-52461.5533707477];
    absError = Beci(:)-BeciSol(:);
    relError = absError./Beci(:);
    relError(~isfinite(relError)) = 0;   % Sets Inf and NaN to 0
    same = any( (abs(absError) > absTol) | (abs(relError) > relTol) );
    if same
        errorCount = errorCount +1;
        fprintf(1,'magIgrf error! Beci not BeciSol\n');
        Beci
        BeciSol
    else
        fprintf(1,'Test %d magIgrf OK!\n',testCount);
    end
    
    % Official WMM 2015 test values
    reciSol = [-4952.5257175308;-4175.9667490388;7.0511621387];
    Beci = magWmm(reciSol,timeSol); % TEST
    testCount = testCount + 1;
    BeciSol = [-7946.3504665208;-7225.2267426330;37546.9281707674];
    absError = Beci(:)-BeciSol(:);
    relError = absError./Beci(:);
    relError(~isfinite(relError)) = 0;   % Sets Inf and NaN to 0
    same = any( (abs(absError) > absTol) | (abs(relError) > relTol) );
    if same
        errorCount = errorCount +1;
        fprintf(1,'magWmm error! Beci not BeciSol\n');
        Beci
        BeciSol
    else
        fprintf(1,'Test %d magWmm OK!\n',testCount);
    end
    
    % Magnetic Field Calculator http://www.ngdc.noaa.gov/geomag-web/?model=igrf#igrfwmm
    Beci = magIgrf(reciSol,timeSol); % TEST
    testCount = testCount + 1;
    BeciSol = [-7948.1795179286;-7222.5871051027;37549.9309835582];
    absError = Beci(:)-BeciSol(:);
    relError = absError./Beci(:);
    relError(~isfinite(relError)) = 0;   % Sets Inf and NaN to 0
    same = any( (abs(absError) > absTol) | (abs(relError) > relTol) );
    if same
        errorCount = errorCount +1;
        fprintf(1,'magIgrf error! Beci not BeciSol\n');
        Beci
        BeciSol
    else
        fprintf(1,'Test %d magIgrf OK!\n',testCount);
    end
    
    % Official WMM 2015 test values
    timeSol.year=2017;
    timeSol.mon = 7;
    timeSol.day = 3;
    timeSol.hr =  0;
    timeSol.min= 0;
    timeSol.sec= 0.0;
    reciSol = [222.3186744526;-1091.0602700267;6259.1315168516];
    Beci = magWmm(reciSol,timeSol); % TEST
    testCount = testCount + 1;
    BeciSol = [-3440.7027887425;15605.1788907721;-52479.4917421705];
    absError = Beci(:)-BeciSol(:);
    relError = absError./Beci(:);
    relError(~isfinite(relError)) = 0;   % Sets Inf and NaN to 0
    same = any( (abs(absError) > absTol) | (abs(relError) > relTol) );
    if same
        errorCount = errorCount +1;
        fprintf(1,'magWmm error! Beci not BeciSol\n');
        Beci
        BeciSol
    else
        fprintf(1,'Test %d magWmm OK!\n',testCount);
    end
    
    % Magnetic Field Calculator http://www.ngdc.noaa.gov/geomag-web/?model=igrf#igrfwmm
    Beci = magIgrf(reciSol,timeSol); % TEST
    testCount = testCount + 1;
    BeciSol = [-3440.1868719256;15606.2282269078;-52475.5709356830];
    absError = Beci(:)-BeciSol(:);
    relError = absError./Beci(:);
    relError(~isfinite(relError)) = 0;   % Sets Inf and NaN to 0
    same = any( (abs(absError) > absTol) | (abs(relError) > relTol) );
    if same
        errorCount = errorCount +1;
        fprintf(1,'magIgrf error! Beci not BeciSol\n');
        Beci
        BeciSol
    else
        fprintf(1,'Test %d magIgrf OK!\n',testCount);
    end
    
    % Official WMM 2015 test values
    reciSol = [4890.0569722041;4248.9453834772;-8.0577661365];
    Beci = magWmm(reciSol,timeSol); % TEST
    testCount = testCount + 1;
    BeciSol = [7900.3409632929;7085.5154738892;37572.5445674497];
    absError = Beci(:)-BeciSol(:);
    relError = absError./Beci(:);
    relError(~isfinite(relError)) = 0;   % Sets Inf and NaN to 0
    same = any( (abs(absError) > absTol) | (abs(relError) > relTol) );
    if same
        errorCount = errorCount +1;
        fprintf(1,'magWmm error! Beci not BeciSol\n');
        Beci
        BeciSol
    else
        fprintf(1,'Test %d magWmm OK!\n',testCount);
    end
    
    % Magnetic Field Calculator http://www.ngdc.noaa.gov/geomag-web/?model=igrf#igrfwmm
    Beci = magIgrf(reciSol,timeSol); % TEST
    testCount = testCount + 1;
    BeciSol = [7906.7937764379;7094.7028926387;37570.0340815931];
    absError = Beci(:)-BeciSol(:);
    relError = absError./Beci(:);
    relError(~isfinite(relError)) = 0;   % Sets Inf and NaN to 0
    same = any( (abs(absError) > absTol) | (abs(relError) > relTol) );
    if same
        errorCount = errorCount +1;
        fprintf(1,'magIgrf error! Beci not BeciSol\n');
        Beci
        BeciSol
    else
        fprintf(1,'Test %d magIgrf OK!\n',testCount);
    end
    
    % My birthday in three years 5km above Trondheim (lat=63, lon=10, alt=5)
	timeSol.year=2018;
	timeSol.mon = 10;
	timeSol.day = 15;
	timeSol.hr =  12;
	timeSol.min= 34;
	timeSol.sec= 56.789;
	reciSol = [-2132.3132036918;-1962.8765465817;5668.2203696986];
	tic;
	Beci = magWmm(reciSol,timeSol); % TEST
    toc;
	testCount = testCount + 1;
	BeciSol = [26150.4541807207;22986.7162207562;-38185.5741861825];
	same = any( (abs(absError) > absTol) | (abs(relError) > relTol) );
    if same
        errorCount = errorCount +1;
        fprintf(1,'magWmm error! Beci not BeciSol\n');
        Beci
        BeciSol
    else
        fprintf(1,'Test %d magWmm OK!\n',testCount);
    end

	tic;
	Beci = magIgrf(reciSol,timeSol); % TEST
	toc;
	testCount = testCount + 1;
	BeciSol = [26166.5929731794;22972.4009119709;-38220.2463841567];
	if same
        errorCount = errorCount +1;
        fprintf(1,'magIgrf error! Beci not BeciSol\n');
        Beci
        BeciSol
    else
        fprintf(1,'Test %d magIgrf OK!\n\n',testCount);
    end
    
%% Test summary
    fprintf(1,'There is(are) %d error(s) in the %d tests.\n',errorCount,testCount);
    