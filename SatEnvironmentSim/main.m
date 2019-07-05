

%%==============================================================================
% Notes
%%==============================================================================
%
%
%
%    
%    Make SPICE a part of this simulation. use it to replace the parts meant to
%    runn on the satelite. Usefull to learn probobaly
%    
%    
%

%% Code to make Octave behave as desired 
addpath( genpath( "D:/Users/UrgST/Documents/GitRepos/test",".git" ) )
format long

%% Code to make Matlab behave as desired
% PLACEHOLDER
%%

%% =============================================================================
%% Sat Orbit environment Sim
%% =============================================================================



clear all

%---------------------------- Declaring Constants ------------------------------
run constants.m
AU = 149597870700; % Astronomical units in [m]
km2m = 1000;

%-----------------------Initial Earth reflection map data ----------------------

refl = resize_refl( load( 'ga050101-051231.mat' ), 2 );
refl.radiMap = idx2RadiMap( refl.data );
refl.areaMap = sphrFaceAreas( refl.data, refl.radiMap );
refl.normalMap = sphrFaceNormals( refl.data );

%---------------------- Initializing satelite orbit ----------------------------

satrec = orbitTwoline2rv(tle_testSat);

%-----------------------------Simulation Time-----------------------------------  3 hours = 56.6MB  24h = 452MB
%                                               year,mon,day,hour,min,sec
[Start_Jd, End_Jd, Delta_Jd] = Init_Jd_stepping(1980, 10, 01, 01, 30, 00.0,...  %Start
                                                1980, 10, 01, 04, 30, 00.0,...  %End
                                                                  00, 01.0);    %steps size

%------------------------Initializing storage struct ---------------------------

steps.Jd = [];
steps.a_sat_ecef = struct( 'irr', 'vect' );
steps.eci2ecef = zeros( 3 );
steps.r_sun_ecef = zeros( 3, 1 );
steps.r_sat_ecef = zeros( 3, 1 );
steps.B_sat_ecef = zeros( 3, 1 );

S = ceil( (End_Jd - Start_Jd) / Delta_Jd );
steps(1:S) = steps;
s = 0;



for Jd = Start_Jd:Delta_Jd:End_Jd;
        
        
        timeUTC = timeJd2datetime( Jd );
        timeUTC.sec = round( timeUTC.sec * 100 ) / 100;
        
        %----------------------Satelite Position Velocity in ecef---------------
        [satrec2,r_sat_eci,~] = orbitSgp4( satrec, timeUTC );
        eci2ecef = frameEci2ecef3( timeUTC, false );
        r_sat_ecef = eci2ecef * r_sat_eci * km2m;
        
        %-------------------Environmental effects position dependent------------
        r_sun_eci = sun( timeUTC ) * AU;
        r_sun_ecef = eci2ecef * r_sun_eci;
        
        a_sat_ecef = albedoV( r_sat_ecef, r_sun_ecef,refl );
        
        B_sat_ecef = magIgrf( r_sat_ecef, timeUTC );
        
        %gravity vector      variations to small to matter?
        
        %aerodynamic torque  size and shape of satelite
        
        %solar prezure       size and shape of satelite
        
        %-----------------------------------------------------------------------
                
         %frame transform to nadir frame frame
                %definition z axis is nadir pointing. y axis is in the negative orbit normal direction(the frame rotates acording to the left hand rule alonge the y axis), the x axis is a result of these two and in circular orbit cases is in orbit velocity direction
                %this is the frame as defined in the nuts project, instead of orbit frame nadir frame seems like a more apropriate name. It seems and orbit frame would more intuitively be one with z axis pointing in the orbit direction x axis in orbit velocity direction and y axis as a result of these
                %the nadir frame would in many cases be usefull for pointing towards earth and orbit frame usefull for spacecraft orbit controll
                
        %----------------------------------- store in --------------------------

        s += 1;
        steps(s).Jd = Jd;
        steps(s).eci2ecef = eci2ecef;
        steps(s).r_sat_ecef = r_sat_ecef;
        steps(s).r_sun_ecef = r_sun_ecef;
        steps(s).a_sat_ecef = a_sat_ecef;
        steps(s).B_sat_ecef = B_sat_ecef;
        
end

% ------------------------------Save simulation --------------------------------

satname = "aTestDay";
starttime = "_1980y10m01d01h30m00s";
Ssteps = num2str(S);
stepsize = ["D10s"];
reflectionmap = "_a050101-051231";
filename = [satname,starttime,Ssteps,stepsize,reflectionmap]

save("-v6", filename,"steps","refl")


%This is where you should place your data illustarator
