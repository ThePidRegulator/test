% IDX2RAD Transform TOMS REFL matrix indices to radians.
%
% [theta,phi] = idx2rad(i,j,sy,sx)
%
% $Id: idx2rad.m,v 1.5 2006/05/17 14:39:18 danji Exp $

function [theta,phi] = idx2rad(i,j,sy,sx);

CONST.d2r = pi/180;

dx = 2*pi/sx;
dy = pi/sy;

phi = pi-dy/2-(i-1)*dy; % y: 1 -> sy, phi: 180 -> 0
theta = (j-1)*dx-pi+dx/2; % x: 1-> sx theta: -180 -> 180

return
