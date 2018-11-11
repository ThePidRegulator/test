%Orbit Position Prediction
% function orbitSgp4
%
% Created: 03.02.2015 15:16:27
% Author: Antoine Pignede
%
% This procedure is the sgp4 prediction model from space command. this is an
%   updated and combined version of sgp4 and sdp4, which were originally
%   published separately in spacetrack report #3. this version follows the
%   methodology from the aiaa paper (2006) describing the history and
%   development of the code.
%
% Copied and edited from Vallado (2013) Fundamentals of Astrodynamics and Applications
%
%   inputs        :
%    satrec      - initialised structure from sgp4init() call.
%    time        - struct containing:
%      year        - year                           1900 .. 2100
%      mon         - month                          1 .. 12
%      day         - day                            1 .. 28,29,30,31
%      hr          - universal time hour            0 .. 23
%      min         - universal time min             0 .. 59
%      sec         - universal time sec             0.0 .. 59.999
%
%   outputs       :
%     reci       - eci position vector                     km
%     veci       - eci velocity vector                     km/s
%
%   locals        :
%     am          -
%     axnl, aynl        -
%     betal       -
%     COSIM   , SINIM   , COSOMM  , SINOMM  , Cnod    , Snod    , Cos2u   ,
%     Sin2u   , Coseo1  , Sineo1  , Cosi    , Sini    , Cosip   , Sinip   ,
%     Cosisq  , Cossu   , Sinsu   , Cosu    , Sinu
%     Delm        -
%     Delomg      -
%     Ecose       -
%     El2         -
%     Eo1         -
%     Esine       -
%     Argpm       -
%     Argpp       -
%     Pl          -
%     Rl          -
%     Su          -
%     T2  , T3   , T4    , Tc
%     Tem5, Temp , Temp1 , Temp2  , Tempa  , Tempe  , Templ
%     U   , Ux   , Uy    , Uz     , Vx     , Vy     , Vz
%     inclm       - inclination
%     mm          - mean anomaly
%     nm          - mean motion
%     nodem      - longi of ascending node
%     xinc        -
%     xincp       -
%     xl          -
%     xlm         -
%     mp          -
%     xmdf        -
%     xmx         -
%     xmy         -
%     nodedf     -
%     xnode       -
%     nodep      -
%
%   coupling      :
%     timeDatetime2jd   - convert date and time values to julian date
%     frameTeme2eci     - convert position in teme from sgp4 to eci
%
% Note that the original function also has all calculations for deep space
%   objects. Those were removed here because this program only will run on
%   a satellite in low earth orbit. The velocity is not needed either.
%   In contrast to the original function, the user has to provide the
%   desired time and date for the position. The time since satellite epoch
%   in minutes is calculated first in this function.
%
% See also
%   orbitTwoline2rv.m, orbitSgp4init.m  - functions to initialize satrec

function [satrec, reci, veci] = orbitSgp4(satrec,time) %#codegen
global TWOPI REQU VELKMPS XKE J2

   % tsince    - time since epoch (minutes)
   jd = timeDatetime2jd(time);
   tsince = (jd - satrec.jdsatepoch) * 1440;
   satrec.t     = tsince;
   
   mrt = 0.0;

   % /* ------- update for secular gravity and atmospheric drag ----- */
   xmdf    = satrec.mo + satrec.mdot * satrec.t;
   argpdf  = satrec.argpo + satrec.argpdot * satrec.t;
   nodedf  = satrec.nodeo + satrec.nodedot * satrec.t;
   argpm   = argpdf;
   mm      = xmdf;
   t2      = satrec.t * satrec.t;
   nodem   = nodedf + satrec.nodecf * t2;
   tempa   = 1.0 - satrec.cc1 * satrec.t;
   tempe   = satrec.bstar * satrec.cc4 * satrec.t;
   templ   = satrec.t2cof * t2;

   if (satrec.isimp ~= 1)
       delomg = satrec.omgcof * satrec.t;
       delm   = satrec.xmcof *...
           ((1.0 + satrec.eta * cos(xmdf))^3 -...
           satrec.delmo);
       temp   = delomg + delm;
       mm     = xmdf + temp;
       argpm  = argpdf - temp;
       t3     = t2 * satrec.t;
       t4     = t3 * satrec.t;
       tempa  = tempa - satrec.d2 * t2 - satrec.d3 * t3 -...
           satrec.d4 * t4;
       tempe  = tempe + satrec.bstar * satrec.cc5 * (sin(mm) -...
           satrec.sinmao);
       templ  = templ + satrec.t3cof * t3 + t4 * (satrec.t4cof +...
           satrec.t * satrec.t5cof);
   end

   nm    = satrec.no;
   em    = satrec.ecco;
   inclm = satrec.inclo;

   if (nm <= 0.0)
       satrec.error = sprintf('mean motion %f less than 0.0',nm);
       reci = [0;0;0];
       veci = [0;0;0];
       return
   end
   
   am = (XKE / nm)^(2.0 / 3.0) * tempa * tempa;
   em = em - tempe;

   % fix tolerance for error recognition
   if ((em >= 1.0) || (em < -0.001) || (am < 0.95))
       satrec.error = sprintf('mean elements, ecc %f >= 1.0 or ecc %f < -0.001 or a %f < 0.95',em,em,am);
       reci = [0;0;0];
       veci = [0;0;0];
       return
   end
   
   % sgp4fix change test condition for eccentricity
   if (em < 1.0e-6)
       em  = 1.0e-6;
   end
   mm     = mm + satrec.no * templ;
   xlm    = mm + argpm + nodem;
   nodem  = rem(nodem, TWOPI);
   argpm  = rem(argpm, TWOPI);
   xlm    = rem(xlm, TWOPI);
   mm     = rem(xlm - argpm - nodem, TWOPI);

   % /* ----------------- compute extra mean quantities ------------- */
   sinim = sin(inclm);
   cosim = cos(inclm);

   % /* -------------------- add lunar-solar periodics -------------- */
   ep     = em;
   xincp  = inclm;
   argpp  = argpm;
   nodep  = nodem;
   mp     = mm;
   sinip  = sinim;
   cosip  = cosim;

   % /* -------------------- long period periodics ------------------ */
   axnl = ep * cos(argpp);
   temp = 1.0 / (am * (1.0 - ep * ep));
   aynl = ep* sin(argpp) + temp * satrec.aycof;
   xl   = mp + argpp + nodep + temp * satrec.xlcof * axnl;

   % /* --------------------- solve kepler's equation --------------- */
   u    = rem(xl - nodep, TWOPI);
   eo1  = u;
   tem5 = 9999.9;
   ktr = 1;
   % sgp4fix for kepler iteration
   % the following iteration needs better limits on corrections
   while (( abs(tem5) >= 1.0e-12) && (ktr <= 10) )
       sineo1 = sin(eo1);
       coseo1 = cos(eo1);
       tem5   = 1.0 - coseo1 * axnl - sineo1 * aynl;
       tem5   = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
       if(abs(tem5) >= 0.95)
           if tem5 > 0.0
               tem5 = 0.95;
           else
               tem5 = -0.95;
           end
       end
       eo1    = eo1 + tem5;
       ktr = ktr + 1;
   end

   % /* ------------- short period preliminary quantities ----------- */
   ecose = axnl*coseo1 + aynl*sineo1;
   esine = axnl*sineo1 - aynl*coseo1;
   el2   = axnl*axnl + aynl*aynl;
   pl    = am*(1.0-el2);
   if (pl < 0.0)
       satrec.error = sprintf('semi-latus rectum %f < 0.0',pl);
       reci = [0;0;0];
       veci = [0;0;0];
       return
   else
       rl     = am * (1.0 - ecose);
       rdotl  = sqrt(am) * esine/rl;
       rvdotl = sqrt(pl) / rl;
       betal  = sqrt(1.0 - el2);
       temp   = esine / (1.0 + betal);
       sinu   = am / rl * (sineo1 - aynl - axnl * temp);
       cosu   = am / rl * (coseo1 - axnl + aynl * temp);
       su     = atan2(sinu, cosu);
       sin2u  = (cosu + cosu) * sinu;
       cos2u  = 1.0 - 2.0 * sinu * sinu;
       temp   = 1.0 / pl;
       temp1  = 0.5 * J2 * temp;
       temp2  = temp1 * temp;

       % /* -------------- update for short period periodics ------------ */
       mrt   = rl * (1.0 - 1.5 * temp2 * betal * satrec.con41) +...
           0.5 * temp1 * satrec.x1mth2 * cos2u;
       su    = su - 0.25 * temp2 * satrec.x7thm1 * sin2u;
       xnode = nodep + 1.5 * temp2 * cosip * sin2u;
       xinc  = xincp + 1.5 * temp2 * cosip * sinip * cos2u;
       mvt   = rdotl - nm * temp1 * satrec.x1mth2 * sin2u / XKE;
       rvdot = rvdotl + nm * temp1 * (satrec.x1mth2 * cos2u +...
           1.5 * satrec.con41) / XKE;

       % /* --------------------- orientation vectors ------------------- */
       sinsu =  sin(su);
       cossu =  cos(su);
       snod  =  sin(xnode);
       cnod  =  cos(xnode);
       sini  =  sin(xinc);
       cosi  =  cos(xinc);
       xmx   = -snod * cosi;
       xmy   =  cnod * cosi;
       ux    =  xmx * sinsu + cnod * cossu;
       uy    =  xmy * sinsu + snod * cossu;
       uz    =  sini * sinsu;
       vx    =  xmx * cossu - cnod * sinsu;
       vy    =  xmy * cossu - snod * sinsu;
       vz    =  sini * cossu;

       % /* -------------------- position (in km) ---------------------- */
       rteme = [0;0;0];
       rteme(1) = (mrt * ux)* REQU;
       rteme(2) = (mrt * uy)* REQU;
       rteme(3) = (mrt * uz)* REQU;
       reci = frameTeme2eci(rteme,time);
       vteme = [0;0;0];
       vteme(1) = (mvt * ux + rvdot * vx) * VELKMPS;
       vteme(2) = (mvt * uy + rvdot * vy) * VELKMPS;
       vteme(3) = (mvt * uz + rvdot * vz) * VELKMPS;
       veci = frameTeme2eci(vteme,time);
   end

   % sgp4fix for decaying satellites
   if (mrt < 1.0)
       satrec.error = sprintf('satellite has decayed, mrt %f < 1.0',mrt);
   end
end
