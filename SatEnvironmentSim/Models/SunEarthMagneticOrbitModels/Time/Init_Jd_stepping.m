

%  
%  inputs        : description                   range / units
%      yr        - year                          1900 .. 2100
%      mn        - month                         1 .. 12
%      d         - day                           1 .. 28,29,30,31
%      h         - universal time hour           0 .. 23
%      m         - universal time min            0 .. 59
%      s         - universal time sec            0.0 .. 59.999
%
%  outputs       :
%    Jd          - julian date                    days from 4713 bc
%
%  Prenominators :
%     B          - Begining
%     E          - End
%     D          - delta

%  locals        :
%    none.
%
%  coupling      :
%    timeDatetime2jd.m


function [B_Jd, E_Jd, D_Jd] = Init_Jd_stepping( B_yr,B_mn,B_d,B_h,B_m,B_s,...
                                                E_yr,E_mn,E_d,E_h,E_m,E_s,...
                                                                  D_m,D_s);

                                                                        
    
    
    %start time struct
    start_timeUTC.year  =   B_yr;
    start_timeUTC.mon   =   B_mn;
    start_timeUTC.day   =    B_d;
    start_timeUTC.hr    =    B_h;
    start_timeUTC.min   =    B_m;
    start_timeUTC.sec   =    B_s;

    %sim stop time struct
    stop_timeUTC.year  =   E_yr;
    stop_timeUTC.mon   =   E_mn;
    stop_timeUTC.day   =    E_d;
    stop_timeUTC.hr    =    E_h;
    stop_timeUTC.min   =    E_m;
    stop_timeUTC.sec   =    E_s;

    %1 second in Jd format
    Jd_sec = 1/86400;                                                             %there is precisely 86400 seconds in a Julian day giving the base unit of 1 second in Jd format

    D_Jd  = Jd_sec*D_s + Jd_sec*60*D_m;                                          
    E_Jd  = timeDatetime2jd(stop_timeUTC);
    B_Jd  = timeDatetime2jd(start_timeUTC);
  
endfunction