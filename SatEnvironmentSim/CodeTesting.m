

%load("aTestDay_1980y10m01d01h30m00s10800D10s_a050101-051231V2")
%
IllustratetOrbitEnvironment([steps(:).a_sat_ecef], {steps(:).r_sat_ecef}, {steps(:).r_sun_ecef}, refl)
%
%tic
%fovMap2 = earthfovV( [steps(:).r_sun_ecef], refl.normalMap );
%toc