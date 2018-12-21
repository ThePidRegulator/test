%What is Main going to contain
%Description of what main dose 
%Main will be where different orbit cenarios will can be defined
%Once a orbit scenario is simulated it will be stored in a format that allows for
%faster SatASDC simulation
%This should be capable of simulating:
%Other environmental torquers on Orientation

%goal number one get the old sun sensor sim upp

%how much should i hide in this main function

%the satelite needs to
%SENS
%magvect
%sunvect
%rotation velocity
%DISTURBANCE
%earth_refl_sunvect
%aerodynamic torque
%gravity
%solar pressure
%SIMULATION FACTS
%timeUTC actual
%
%

addpath(genpath("D:/Users/UrgST/Documents/GitRepos/test",".git"))
run constants.m

%initiating sim start time struct
start_timeUTC.year  =   1980;
start_timeUTC.mon   =      4;
start_timeUTC.day   =      4;
start_timeUTC.hr    =      1;
start_timeUTC.min   =      0;
start_timeUTC.sec   =    24.0000;

%initiating sim stop time struct
stop_timeUTC.year  =   1980;
stop_timeUTC.mon   =      4;
stop_timeUTC.day   =      20;
stop_timeUTC.hr    =      2;
stop_timeUTC.min   =      1;
stop_timeUTC.sec   =   24.00;


%initiating step size is d_jd = step_time_Jd-start_time_Jd
step_timeUTC.year  =      1980;
step_timeUTC.mon   =      4;
step_timeUTC.day   =      4;
step_timeUTC.hr    =      1;
step_timeUTC.min   =      20;
step_timeUTC.sec   =      24.00;

d_Jd           = timeDatetime2jd(step_timeUTC)-timeDatetime2jd(start_timeUTC);      % + 0.00000627779999e-05; %this stabelizes the steping in seconds for many years
stop_time_Jd   = timeDatetime2jd(stop_timeUTC);
start_time_Jd  = timeDatetime2jd(start_timeUTC);

satrec = orbitTwoline2rv(tle);

for Jd = start_time_Jd:d_Jd:stop_time_Jd

  timeUTC                       = timeJd2datetime(Jd);                              % may need a stabelizing function for given step d_Jd, conversion of 1 second to Jd is not precise enouf for stepping;
  timeUTC.sec                   = round(timeUTC.sec*10)/10;                         % this stabelizes the stepping 
  Jd                            = timeDatetime2jd(timeUTC);
  
  [satrec2,r_sat_eci,v_sat_eci]   = orbitSgp4(satrec,timeUTC);
  
  Beci                          = magIgrf(r_sat_eci,timeUTC);
  
  r_sun_eci                     = sun(timeUTC);
  r_sun_ecef                    = frameEci2ecef(r_sun_eci,timeUTC);
  r_sat_ecef                    = frameEci2ecef(r_sat_eci,timeUTC);
  
  %albedo                        = albedo_wrapper(r_sat_ecef,r_sun_ecef,param,redfac,refllib);
  %albedo_v                      = CreateAlbedoVectors(r_sat_ecef,Albedo_ecef);
  
  %gravity vector      variations to small to matter
  
  %aerodynamic torque  size and shape of satelite
  
  %solar prezure       size and shape of satelite
  
  timeUTC.mon
  timeUTC.day
  timeUTC.min
  
end


