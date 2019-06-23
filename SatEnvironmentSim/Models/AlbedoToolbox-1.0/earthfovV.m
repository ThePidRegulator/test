% A = earthfovV2()
% 
% Created: 21.04.2019 16:35:00
% Author: Henrik Rudi Haave
%
% Finds faces of a equirectangular maped Earth/sphere that can be seen from a
% or more positions above, and returned as logical maps. For one position the 
% function simply tests all faces for visibility using the folowing equation.
%
%        faceNormal * (pos - facePos) >= 0
% faceNormal * (pos - faceNormal*EMR) >= 0, faceNormal*faceNormal = 1
%                      faceNormal*pos >= EMR
%
% For many positions this equation is only used on a rectangular subset of 
% faces. This function is faster than the old earthfov that was based on 
% trigonometric functions.
%
%  Inputs        : description                                 format / units
%             pos- N positions above Earth                     [](3xN)/ meters
%    facesNormMap- normal vectors of faces SY < SX         [](SYxSXx3)/ 1 unit
%                -  
%  Outputs       :
%         fovMaps- maps of viewable areas from pos         [](SYxSXxN)/ logical
%  Locals        :
%       CONST.EMR- mean radius or Earth                               / meters
%    facesNormMap- normal vectors of sphere/Earth faces               / UnitVect
%            CY2-  radian map increment size N..S                     / radians
%            CX2-  radian map increment size W..S                     / radians
%             SY-  N..S total map elements/bins/areas                 / +integer
%             SX-  W..E total map elements/bins/areas                 / +integer
%
%  Coupling      :
%
%  See also      :
%

function fovMaps = earthfovV( pos, facesNormMap );
        CONST.EMR = 6371.01e3;
        
        N = size(pos,2);
        %Make logic map of visible sphare faces for 1 position
        if N == 1
                pos = reshape( pos, [1,1,3] );
                fovMaps = sum( facesNormMap .* pos , 3 ) >= CONST.EMR;
                return
        end
        
        % Find sub rectangles of map for reduced amount of faces to test 
        [SY, SX, ~] = size( facesNormMap );
        CX2 = 2*pi / SX;
        CY2 = pi / SY;
    
        posSph = cart2sph( pos(1,:),pos(2,:),pos(3,:) );
        posSph(:, 2) = pi/2 - posSph(:, 2);%phi, spher to polar
        rho =  acos( CONST.EMR ./ posSph(:, 3) );
        
        idxOffsetX = ceil( rho./(sin( posSph(:, 2) ).*CX2) );
        idxOffsetY = ceil( rho./(CY2) );
        
        idxPosX = round( (posSph(:, 1) + pi) ./ CX2 );
        idxPosY = round( (pi-posSph(:, 2)) ./ CY2 );
  
        fovMinX = idxPosX - idxOffsetX;
        fovMaxX = idxPosX + idxOffsetX;
        
        fovMinY = max( 1, idxPosY - idxOffsetY );
        fovMaxY = min( SY, idxPosY + idxOffsetY );

        %Wrapping fovMaxX indexes above SX to above or equal 1 
        wrapMaxX = zeros(N,1);
        wrapMaxX(fovMaxX>SX) = mod( fovMaxX(fovMaxX>SX), SX );
        fovMaxX(fovMaxX>SX) = SX;

        %Wrapping fovMinX indexes belov 1 to belov or equal SX 
        wrapMinX = ones(N,1)*(SX + 1);
        wrapMinX(fovMinX<1) = SX + mod( 1, fovMinX(fovMinX<1)-1 );
        fovMinX(fovMinX<1) = 1;
        
        fovMaxX = num2cell(fovMaxX);
        fovMinX = num2cell(fovMinX);
        fovMinY = num2cell(fovMinY);
        fovMaxY = num2cell(fovMaxY);
        wrapMinX = num2cell(wrapMinX);
        wrapMaxX = num2cell(wrapMaxX);
        
        fovMap = zeros( SY, SX, 'logical' );
        pos = reshape( pos', [N,1,3] );
        pos = num2cell( pos, 3 );
        
        % finding faces/map elements within fov
        fovMaps = cellfun( @posfov...
                , {facesNormMap}...
                , {fovMap}...
                , {CONST.EMR}...
                , pos...
                , fovMinY, fovMaxY...
                , wrapMaxX, fovMinX, fovMaxX, wrapMinX, {SX}...
                , 'UniformOutput', false ...
        );
        % converting from cell to matrix, this is slow and takes 0.3 sec for 10800 maps
        fovMaps = cat( 3, fovMaps{:} );
        
end

function fovMap = posfov( facesNormMap, fovMap,EMR, pos, x,y,a,z,w,b,SX )
      
        X = [1:a z:w b:SX];
        Y = x:y;
        
        fovMap(Y,X) = sum( facesNormMap(Y,X,:) .* pos, 3 ) >= EMR;
        
end

% Notes on atempts at making this function fast in Octave
%     Without cellfun the fovMaps equation took 1 sec for 10800 positions 
%     with the code belov, linIdx are indexes of subrectangles:
%     fovMaps(linIdx(:,1)) = sum(facesNormMap(linIdx).*posV , 2) >= EMR 
%     cellfun used as in this function takes aprox 1.65 sec. 
%     concatenation from cell to matrice takes an aditional 0.3 sec.
%     A method for making linIdx and posV faster than a second was not found
%     and used in total atleast 1.5 seconds. Because the "colon" operator can 
%     not be vectorised it is likely imposible to make linIdx without cellfun
%     or a for loop making the proces slow compared to a fully vectorised 
%     operation
%

