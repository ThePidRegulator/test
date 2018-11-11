%Orbit Position Prediction
% function orbitTwoline2rv
%
% Created: 03.02.2015 11:21:25
% Author: Antoine Pignede
%
% This function converts the two line element set character string data to
%   variables and initializes the sgp4 variables. several intermediate varaibles
%   and quantities are determined. simple checks of the input tle are done.
%
% Copied and edited from Vallado (2013) Fundamentals of Astrodynamics and Applications
%
% inputs        :
%   tle         - TLE string in one line without any kind of space between the first and second line
% 
% outputs       :
%   satrec      - structure containing all the sgp4 satellite information
%     The first satrec struct element is an error message. If everything is
%     fine it will be empty, otherwise a short description of the error is
%     given. In orbitTwoline2rv possible errors in the tle string are
%     detected. In orbitSgp4 values outside the expected range are
%     detected. However the calculations will always continue. The only
%     error to stop a function is if the TLE input string in
%     orbitTwoline2rv is less than 138 characters long.
% 
% coupling      :
%   timeDays2datetime  - conversion of days to month, day, hour, minute, second
%   timeDatetime2jd - convert day month year hour minute second to julian date
%   orbitSgp4init   - initialize the sgp4 variables
%
% Note that the original function takes the two-line element set (TLE) as
%   two strings representing one line each. Here the input is just one
%   string containing the TLE without any kind of separation between the
%   two lines. As they always have 69 characters, the string of length 138
%   is separated in the middle.
% The original function uses a special extended TLE to also have the time
%   points where values are to be calculated. Anything related to that was
%   removed here. This function just converts the "normal" TLE to values
%   for the satellite and calls orbitSgp4init to calculate more values.
% This function is not meant to be used to track several satellites.
%   Some easy checks for right TLE input were added.
% 
% See also
%   orbitSgp4.m  - calculate orbital position of satrec object

function [satrec] = orbitTwoline2rv(tle) %#codegen
global DEG2RAD
    
    satrec.error = '';

    % check for right length of input string
    if length(tle) < 138
        satrec.error = sprintf('too short TLE, length(tle) = %d < 138',length(tle));
        return;
    end
    
    % split input string into the two lines
    longstr1 = tle(1:69);
    longstr2 = tle(70:end);

    % set the implied decimal points since doing a formated read 
    % fixes for bad input data values (missing, ...)
    if (longstr1(8) == ' ')
        longstr1(8) = 'U';
    end
    if (longstr1(10) == ' ')
        longstr1(10) = '.';
    end
    for j = 11:16
        if (longstr1(j) == ' ')
            longstr1(j) = '_';
        end
    end
    if (longstr1(45) ~= ' ')
        longstr1(44) = longstr1(45);
    end
    longstr1(45) = '.';     
    for j = 46:50
        if (longstr1(j) == ' ')
            longstr1(j) = '0';
        end
    end    
    if (longstr1(52) == ' ')
        longstr1(52) = '0';
    end    
    if (longstr1(54) ~= ' ')
        longstr1(53) = longstr1(54);
    end    
    longstr1(54) = '.';
    if (longstr1(63) == ' ')
        longstr1(63) = '0';
    end
    if ((length(longstr1) < 68) || (longstr1(68) == ' '))
        longstr1(68) = '0';
    end

    longstr2(26) = '.';     
    for j = 27:33
        if (longstr2(j) == ' ')
            longstr2(j) = '0';
        end
    end

    % parse first line
    linenum1 = str2double(longstr1(1));
    satrec.satnum = str2double(longstr1(3:7));
    satrec.classification = longstr1(8);
    satrec.intldesg = longstr1(10:17);
    satrec.epochyr = str2double(longstr1(19:20));
    satrec.epochdays = str2double(longstr1(21:32));
    satrec.ndot = str2double(longstr1(34:43));
    satrec.nddot = str2double(longstr1(44:50));
    nexp = str2double(longstr1(51:52));
    satrec.bstar = str2double(longstr1(53:59));
    ibexp = str2double(longstr1(60:61));
    satrec.ephemtype = str2double(longstr1(63));
    satrec.elnum = str2double(longstr1(65:68));
    checksum1 = str2double(longstr1(69));
 
    % parse second line
    linenum2 = str2double(longstr2(1));
    satnum2 = str2double(longstr2(3:7));
    satrec.inclo = str2double(longstr2(8:16));
    satrec.nodeo = str2double(longstr2(17:25));
    satrec.ecco = str2double(longstr2(26:33));
    satrec.argpo = str2double(longstr2(34:42));
    satrec.mo = str2double(longstr2(43:51));
    satrec.no = str2double(longstr2(52:63));
    satrec.revnum = str2double(longstr2(64:68));
    checksum2 = str2double(longstr2(69));
    
    % check for line numbers
    if linenum1 ~= 1
        satrec.error = sprintf('wrong first line number, %d not 1',linenum1);
    end
    if linenum2 ~= 2
        satrec.error = sprintf('wrong second line number, %d not 2',linenum2);
    end
    % check for satellite number
    if satrec.satnum ~= satnum2
        satrec.error = sprintf('two different satellite numbers, %d not %d',satrec.satnum,satnum2);
    end
    % check for classification character
    if satrec.classification ~= 'U'
        satrec.error = sprintf('unknown classification, %c not U',satrec.classification);
    end
    % check the first line checksum
    checksum1Calc = 0;
    for i=1:(length(longstr1)-1)
        if longstr1(i) == '-'
            checksum1Calc = checksum1Calc + 1;
        elseif ~isnan(str2double(longstr1(i)))
            checksum1Calc = checksum1Calc + str2double(longstr1(i));
        end
    end
    checksum1Calc = rem(checksum1Calc,10);
    if checksum1Calc ~= checksum1
        satrec.error = sprintf('wrong first line checksum, %d not %d',checksum1Calc,checksum1);
    end
    % check the second line checksum
    checksum2Calc = 0;
    for i=1:(length(longstr2)-1)
        if longstr2(i) == '-'
            checksum2Calc = checksum2Calc + 1;
        elseif ~isnan(str2double(longstr2(i)))
            checksum2Calc = checksum2Calc + str2double(longstr2(i));
        end
    end
    checksum2Calc = rem(checksum2Calc,10);
    if checksum2Calc ~= checksum2
        satrec.error = sprintf('wrong second line checksum, %d not %d',checksum2Calc,checksum2);
    end

    satrec.nddot= satrec.nddot * 10.0^nexp;
    satrec.bstar= satrec.bstar * 10.0^ibexp;
    
    % [rev/day]/[rad/min] = 229.1831180523293;
    xpdotp   =  1440.0 / (2.0*pi);
    
    % ---- find no, ndot, nddot ----
    satrec.no   = satrec.no / xpdotp;                      % rad/min
    satrec.ndot = satrec.ndot  / (xpdotp*1440.0);          % rad/min^2
    satrec.nddot= satrec.nddot / (xpdotp*1440.0*1440);     % rad/min^3

    % ---- find standard orbital elements ----
    satrec.inclo = satrec.inclo  * DEG2RAD;
    satrec.nodeo = satrec.nodeo * DEG2RAD;
    satrec.argpo = satrec.argpo  * DEG2RAD;
    satrec.mo    = satrec.mo     *DEG2RAD;

    % ------------- temp fix for years from 1957-2056 ----------------
    % ------ correct fix will occur when year is 4-digit in 2le ------
     if (satrec.epochyr < 57)
         year= satrec.epochyr + 2000;
       else
         year= satrec.epochyr + 1900;
     end;
     time = timeDays2datetime(year,satrec.epochdays);
     satrec.jdsatepoch = timeDatetime2jd(time);
     
     % ------------- initialize the orbit at sgp4epoch --------------
     [satrec] = orbitSgp4init(satrec);
end
     