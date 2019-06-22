%Orbit Position Prediction
% function orbitSgp4init
%
% Created: 03.02.2015 12:02:50
% Author: Antoine Pignede
%
% This procedure initializes variables for sgp4. It is automatically called
%   each time a new tle has been read.
%
% Copied and edited from Vallado (2013) Fundamentals of Astrodynamics and Applications
%
%   inputs        :
%     satrec      - structure containing all the sgp4 satellite information
%
%   outputs       :
%     satrec      - structure containing all the sgp4 satellite information
%
%   locals        :
%     Cc1sq  , Cc2    , Cc3
%     Coef   , Coef1
%     cosio4      -
%     eeta        -
%     etasq       -
%     perige      - perigee
%     pinvsq      -
%     psisq       -
%     qzms24      -
%     sfour       -
%     temp        -
%     temp1, temp2, temp3       -
%     tsi         -
%     xhdot1      -
%
%   coupling      :
%    none
%
% Note that the original function also has all calculations for deep space
%   objects. Those were removed here because this program only will run on
%   a satellite in low earth orbit. The velocity is not needed either.
%   This function is not used directly but as subroutine of orbitTwoline2rv.m.
%
% See also
%   orbitTwoline2rv.m  - convert TLE string to satrec values

function [satrec] = orbitSgp4init(satrec) %#codegen
global REQU XKE J2 J4 J3OJ2

%% Initialization
   % /* ----------- set all near earth variables to zero ------------ */
   satrec.isimp   = 0;   satrec.aycof  = 0.0;
   satrec.con41   = 0.0; satrec.cc1    = 0.0; satrec.cc4      = 0.0;
   satrec.cc5     = 0.0; satrec.d2     = 0.0; satrec.d3       = 0.0;
   satrec.d4      = 0.0; satrec.delmo  = 0.0; satrec.eta      = 0.0;
   satrec.argpdot = 0.0; satrec.omgcof = 0.0; satrec.sinmao   = 0.0;
   satrec.t       = 0.0; satrec.t2cof  = 0.0; satrec.t3cof    = 0.0;
   satrec.t4cof   = 0.0; satrec.t5cof  = 0.0; satrec.x1mth2   = 0.0;
   satrec.x7thm1  = 0.0; satrec.mdot   = 0.0; satrec.nodedot = 0.0;
   satrec.xlcof   = 0.0; satrec.xmcof  = 0.0; satrec.nodecf  = 0.0;

   ss     = 78.0 / REQU + 1.0;
   qzms2t = ((120.0 - 78.0) / REQU)^4;
   % sgp4fix divisor for divide by zero check on inclination
   % the old check used 1.0 + cos(pi-1.0e-9), but then compared it to
   % 1.5 e-12, so the threshold was changed to 1.5e-12 for consistancy
   temp4    =   1.5e-12;

   satrec.t    = 0.0;

%% Excerpt of initl.m from Vallado (2013) Fundamentals of Astrodynamics and Applications
%  This procedure initializes the spg4 propagator. all the initialization is
%    consolidated here instead of having multiple loops inside other routines.

   % /* ------------- calculate auxillary epoch quantities ---------- */
   eccsq  = satrec.ecco * satrec.ecco;
   omeosq = 1.0 - eccsq;
   rteosq = sqrt(omeosq);
   cosio  = cos(satrec.inclo);
   cosio2 = cosio * cosio;

   % /* ------------------ un-kozai the mean motion ----------------- */
   ak    = (XKE / satrec.no)^(2.0 / 3.0);
   d1    = 0.75 * J2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq);
   del   = d1 / (ak * ak);
   adel  = ak * (1.0 - del * del - del *...
       (1.0 / 3.0 + 134.0 * del * del / 81.0));
   del   = d1/(adel * adel);
   satrec.no    = satrec.no / (1.0 + del);

   ao    = (XKE / satrec.no)^(2.0 / 3.0);
   sinio = sin(satrec.inclo);
   po    = ao * omeosq;
   con42 = 1.0 - 5.0 * cosio2;
   satrec.con41 = -con42-cosio2-cosio2;
   posq  = po * po;
   rp    = ao * (1.0 - satrec.ecco);
   
%% Continue initialization for near space

   if ((omeosq >= 0.0 ) || ( satrec.no >= 0.0))
       satrec.isimp = 0;
       if (rp < (220.0 / REQU + 1.0))
           satrec.isimp = 1;
       end
       sfour  = ss;
       qzms24 = qzms2t;
       perige = (rp - 1.0) * REQU;

       % /* - for perigees below 156 km, s and qoms2t are altered - */
       if (perige < 156.0)
           sfour = perige - 78.0;
           if (perige < 98.0)
               sfour = 20.0;
           end
           qzms24 = ((120.0 - sfour) / REQU)^4.0;
           sfour  = sfour / REQU + 1.0;
       end
       pinvsq = 1.0 / posq;

       tsi  = 1.0 / (ao - sfour);
       satrec.eta  = ao * satrec.ecco * tsi;
       etasq = satrec.eta * satrec.eta;
       eeta  = satrec.ecco * satrec.eta;
       psisq = abs(1.0 - etasq);
       coef  = qzms24 * tsi^4.0;
       coef1 = coef / psisq^3.5;
       cc2   = coef1 * satrec.no * (ao * (1.0 + 1.5 * etasq + eeta *...
           (4.0 + etasq)) + 0.375 * J2 * tsi / psisq * satrec.con41 *...
           (8.0 + 3.0 * etasq * (8.0 + etasq)));
       satrec.cc1   = satrec.bstar * cc2;
       cc3   = 0.0;
       if (satrec.ecco > 1.0e-4)
           cc3 = -2.0 * coef * tsi * J3OJ2 * satrec.no * sinio / satrec.ecco;
       end
       satrec.x1mth2 = 1.0 - cosio2;
       satrec.cc4    = 2.0* satrec.no * coef1 * ao * omeosq *...
           (satrec.eta * (2.0 + 0.5 * etasq) + satrec.ecco *...
           (0.5 + 2.0 * etasq) - J2 * tsi / (ao * psisq) *...
           (-3.0 * satrec.con41 * (1.0 - 2.0 * eeta + etasq *...
           (1.5 - 0.5 * eeta)) + 0.75 * satrec.x1mth2 *...
           (2.0 * etasq - eeta * (1.0 + etasq)) * cos(2.0 * satrec.argpo)));
       satrec.cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 *...
           (etasq + eeta) + eeta * etasq);
       cosio4 = cosio2 * cosio2;
       temp1  = 1.5 * J2 * pinvsq * satrec.no;
       temp2  = 0.5 * temp1 * J2 * pinvsq;
       temp3  = -0.46875 * J4 * pinvsq * pinvsq * satrec.no;
       satrec.mdot     = satrec.no + 0.5 * temp1 * rteosq * satrec.con41 +...
           0.0625 * temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
       satrec.argpdot  = -0.5 * temp1 * con42 + 0.0625 * temp2 *...
           (7.0 - 114.0 * cosio2 + 395.0 * cosio4) +...
           temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
       xhdot1            = -temp1 * cosio;
       satrec.nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) +...
           2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;
       satrec.omgcof   = satrec.bstar * cc3 * cos(satrec.argpo);
       satrec.xmcof    = 0.0;
       if (satrec.ecco > 1.0e-4)
           satrec.xmcof = -(2.0 / 3.0) * coef * satrec.bstar / eeta;
       end
       satrec.nodecf = 3.5 * omeosq * xhdot1 * satrec.cc1;
       satrec.t2cof   = 1.5 * satrec.cc1;

       % sgp4fix for divide by zero with xinco = 180 deg
       if (abs(cosio+1.0) > 1.5e-12)
          satrec.xlcof   = -0.25 * J3OJ2 * sinio *...
              (3.0 + 5.0 * cosio) / (1.0 + cosio);
       else
         satrec.xlcof   = -0.25 * J3OJ2 * sinio *...
              (3.0 + 5.0 * cosio) / temp4;
       end   
       satrec.aycof   = -0.5 * J3OJ2 * sinio;
       satrec.delmo   = (1.0 + satrec.eta * cos(satrec.mo))^3;
       satrec.sinmao  = sin(satrec.mo);
       satrec.x7thm1  = 7.0 * cosio2 - 1.0;

       % /* ----------- set variables if not deep space ----------- */
       if (satrec.isimp ~= 1)
           cc1sq          = satrec.cc1 * satrec.cc1;
           satrec.d2    = 4.0 * ao * tsi * cc1sq;
           temp           = satrec.d2 * tsi * satrec.cc1 / 3.0;
           satrec.d3    = (17.0 * ao + sfour) * temp;
           satrec.d4    = 0.5 * temp * ao * tsi *...
               (221.0 * ao + 31.0 * sfour) * satrec.cc1;
           satrec.t3cof = satrec.d2 + 2.0 * cc1sq;
           satrec.t4cof = 0.25 * (3.0 * satrec.d3 + satrec.cc1 *...
               (12.0 * satrec.d2 + 10.0 * cc1sq));
           satrec.t5cof = 0.2 * (3.0 * satrec.d4 +...
               12.0 * satrec.cc1 * satrec.d3 +...
               6.0 * satrec.d2 * satrec.d2 +...
               15.0 * cc1sq * (2.0 * satrec.d2 + cc1sq));
       end
   end
end
