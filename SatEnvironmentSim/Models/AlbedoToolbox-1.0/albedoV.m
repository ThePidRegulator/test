% ALBEDO Calculation of albedo for a given satellite, Earth, Sun constellation
% and reflectivity data.
%
% a = albedoV(sat, sun, refl [,type])
%
% sat and sun are the vectors to the satelite and Sun in Cartesian Earth Centered 
% Earth Fixed coordinates, respectively. refl is the reflectivity data to use 
% for albedo calculation. Type is 'p' for 3D plots and 's' for spherical. If
% unspecified no plots are generated.
%
% Original verision made by Dan D.V. Bhanderi has been uppdated to use matrix 
% multiplication instead of for loops and trigonometric functions 
% to speed up computation time in Octave
%
%
% 2019 Henrik Rudi Haave

function a = albedo( sat, sun, refl, type )
        
        CONST.EMR = 6371.01e3;
        CONST.AM0 = 1366.9;
        CONST.d2r = pi / 180;
        
        % Support row vectors
        if length(sat) > size(sat,1)
                sat = sat';
        end
        
        if length(sun) > size(sun, 1)
                sun = sun';
        end
        
        % Check for common error, satellite altitude to low 
        if norm(sat) < CONST.EMR
                error('albedo.m: The satellite has crashed into Earth!');
        end
        
        % Data size
        [sy,sx] = size (refl.data);
        
        % Visible elements & Sunlit elements
        visible = earthfovV( sat, refl.normalMap );
        sunlit = earthfovV( sun, refl.normalMap );
        visLit = visible & sunlit;
        
        % sy row
        % sx col
        % Find matrix indexes of sun reflecting Earth/sphr faces seen by sat
        [rowIdx,colIdx] = find( visLit );
        if( (size(rowIdx, 1)>0) && (size(colIdx)>0) )
                
                visLitIdx = sub2ind( [sy, sx], rowIdx, colIdx );
                
                faceNorm(:,1) = refl.normalMap(:,:,1)(visLitIdx);
                faceNorm(:,2) = refl.normalMap(:,:,2)(visLitIdx);
                faceNorm(:,3) = refl.normalMap(:,:,3)(visLitIdx);
                
                facePos = faceNorm * CONST.EMR;
                
                face2sat = sat' - facePos;
                face2satDist = sqrt( dot( face2sat, face2sat, 2 ) );
                face2satDir = rdivide( face2sat, face2satDist );
                
                % Incident irradiance from the sun onto visible faces
                E_in = CONST.AM0 .* refl.areaMap(rowIdx)'...
                                 .* (faceNorm*(sun/norm(sun)));
                
                % Reflected irradiance from earth as seen by the satelite
                a.irr = dot(face2satDir,faceNorm,2) .* ...
                        (E_in .* refl.data(visLitIdx) ./ ...
                        (pi .* face2satDist.^2));
                
                % Reflected irradiance direction
                a.vect = face2satDir;
                
        else
                
                a.vect = [0,0,0];
                a.irr = [0];
                
        end
        
        %-----------------------------------------------------------------------
        % plot the sat visible, sun lit and visLit of these areas on Earth and 
        % plot reflected irradiance as seen by the satelite 
        if nargin > 3
                figure (1);
                subplot (3, 1, 1);
                plot_refl (mask( refl.data, visible ), type, 'no colorbar');
                title( 'Satellite Field of View' );
        end
        
        if nargin > 3
                subplot (3, 1, 2);
                plot_refl( mask( refl.data, sunlit ), type, 'no colorbar' );
                title ( 'Solar Field of View' );
        end
        
        if nargin > 3
                subplot ( 3, 1, 3 );
                plot_refl ( mask( refl.data, visLit ), type, 'no colorbar' );
                title( 'Sunlit Satellite Field of View' );
        end
        
        if nargin > 3
                irradianceMap = zeros (sy, sx);
                irradianceMap(visLit) = a.irr;
                
                figure (2);
                plot_alb (irradianceMap, type);
        end
        %-----------------------------------------------------------------------
return
