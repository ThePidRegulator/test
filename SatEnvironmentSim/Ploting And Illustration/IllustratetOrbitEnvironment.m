
% quiver3 is exrutiatingly slow and consists of 3 line elements. vectors are
% therefor made of Two lines, one short and one long in diferent colors or style to indicate direction

function IllustratetOrbitEnvironment( albedo, rsatEcef, rsunEcef, refl )
	
        
        EMR = 6371.010e3;
        AM0 = 1366.9;
        D2R = pi/180;
        AU2M = 149597870700;
        
        illuData = IllustrationData;
        % sy row
        % sx col
        %% Pass the data neded for illustration to Class illustartionData
        mapSize = size( refl.data );
        illuData.reflMap = refl.data;
        
        % max number of steps, satelite and sun position
        illuData.S = size( rsatEcef, 2 );
        illuData.rSat = rsatEcef;
        illuData.rSun = rsunEcef;
        
        % Albedo directions, strengths and max reflected sunlight
        illuData.albedoDir = {albedo(:).vect};
        illuData.albedoIrr = {albedo(:).irr};
        
        % Prepare products for illustartion
        illuData.satFovMap = earthfovV( [rsatEcef{:}], refl.normalMap );
        illuData.sunlitMap = earthfovV( [rsunEcef{:}], refl.normalMap );
        illuData.litFovMap = illuData.sunlitMap & illuData.satFovMap;
        clear illuData.satFovMap;
        clear illuData.sunlitMap;
        
        % Position of Albedo sources, center of relevant Earth faces
        earthFaceCenters = {refl.normalMap * EMR};
        illuData.albedoSrce = cellfun( @(x,y) {reshape( x([y,y,y]) ,[],3 )} ,...
                earthFaceCenters,...
                num2cell( illuData.litFovMap, [1,2] )...
        );
        
        illuData.albedoMaxIrr = cellfun( @max, ...
                        {albedo(:).irr}, ...
                        'UniformOutput', false ...
        );
        
        
        %% Figure settings
        figure 2, hold on, grid on
        axis("square");
        set( gca, ...
                'Xlim', [-8000000,8000000], ...
                'Ylim', [-8000000,8000000], ...
                'Zlim', [-8000000,8000000], ...
                'color', [0,0,0,], ...
                'gridcolor', [1,1,1]
        );
        
        [sphx,sphy,sphz] = spherMesh( mapSize(1)+1, mapSize(2)+1 );
        EarthH = surf( sphx'*EMR, sphy'*EMR, sphz'*EMR );
        set( EarthH,'Cdata', refl.data' );
        shading flat;
        line( [0;EMR*2], [0;0], [0;0], 'color', 'r', 'linewidth', 3 );
        line( [0;0], [0;EMR*2], [0;0], 'color', 'g', 'linewidth', 3 );
        line( [0;0], [0;0], [0;EMR*2], 'color', 'b', 'linewidth', 3 ); 
        
        line( [0; 3152560]*1.5, [0; 598984]*1.5, [0;5500855]*1.5, 'color', 'white', 'linewidth', 3 ); %Oslo norway
        line( [0; 2302096.036]*1.5, [0; -2223108.300]*1.5, [0; 5505994.590]*1.5, 'color', 'white', 'linewidth', 3 ); %greenland bottom tip
        
        
        illuData.EarthH = EarthH;
        
        % Illustration controll
        sliderH = uicontrol (gcf, ...
                'style', 'slider', ...
                'Units', 'normalized', ...
                'position', [0.001, 0.001, 0.8, 0.05], ...
                'min', 1, ...
                'max', illuData.S, ...
                'callback', {@sliderStepping,illuData},...
                'SliderStep', [1/(illuData.S-1) , 1/(illuData.S-1) ]...
        );
        
        
%       clear "illuSatOrbit"
%       clear "illuSunIrr"
%       clear "illuAlbedoIrr"
%       clf
        %camlookat([rsunH,rsatH])
end



function sliderStepping(hslider ,~ , illuData)
        
        % Get step selected by slider
        s = ceil(get( hslider,"value" ))
        
        illuSteppTo(s, illuData)
        
        drawnow

end

function illuSteppTo(s, illuData)
        
        
        EMR = 6371.01e3;
        AM0 = 1366.9;
        D2R = pi/180;
        AU2M = 149597870700;
        
        rSat = illuData.rSat{s};
        
        scaling = norm(rSat) - EMR;
        
        illuSatOrbit(s, illuData.rSat, 500);
     
        illuSunIrr(rSat, illuData.rSun{s}, scaling);
        
        albedoIrr = illuData.albedoIrr{s};
        maxIrr = illuData.albedoMaxIrr{s};
        albedoDir = illuData.albedoDir{s} * scaling;
        albedoSrce = illuData.albedoSrce{s};
        
        illuAlbedoIrr(albedoSrce, albedoDir, albedoIrr, maxIrr );
        

        reflMapMask = mask( illuData.reflMap, illuData.litFovMap(:,:,s), 0.3 );
        set( illuData.EarthH, 'Cdata', reflMapMask' );
       
end


function illuSatOrbit(s, rSat, sTail)
        
        
        persistent h;
        if isempty(h) || !ishandle(h)
                h = line( rSat{1}(1), rSat{1}(2), rSat{1}(3) );
                set( h,...
                        'linewidth', 2,...
                        'color', [0,1,0]...
                );
        else
                sTail=s-sTail;
                if( sTail<1 )
                        sTail=1;
                end
                
                rSat = reshape( [rSat{sTail:s}], 3, [] );
                
                set( h,...
                        "xdata", rSat(1, :),...
                        "ydata", rSat(2, :),...
                        "zdata", rSat(3, :)...
                );
        end
end
        


function illuSunIrr( position, rSun, scaling) 
        AU2M = 149597870700;
        AM0 = 1366.9;
        
        sun = AM0 * 2000 * -rSun' / AU2M; % rSun points to the sun, reverse it to make it point to coordinate system center, Earth
        sunxyz = position' - sun; % Change position to 'position'
        
        persistent h
        if isempty(h) || !ishandle(h)
                disp('newsunquiver')
                [htipp, hbutt] = fastQuiver(h, sunxyz, sun, 0.2);
                
                set( htipp,...
                        'linewidth', 1.5,...
                        'color', [1,1,0]
                );
                
                set( hbutt,...
                        'linewidth', 1.5,...
                        'color', [1,0.6,0]
                );
                h(1) = htipp;
                h(2) = hbutt;
        else
                fastQuiver(h, sunxyz, sun, 0.2);
        end
        
        
end

function illuAlbedoIrr( albSrc, albDir, albIrr, maxIrr )
        
        maxAlbIrr = max(albIrr);
        if maxAlbIrr > 0; 
                albIrr = albIrr/maxAlbIrr;
        end
        
        persistent h
        if isempty(h) || !ishandle(h)
                [htipp, hbutt] = fastQuiver(h, albSrc, albDir, albIrr);
                set( htipp,...
                        'linewidth', 1,...
                        'color', [1,0,0]
                );
                h(1) = htipp;
                h(2) = hbutt;
        else
                h = fastQuiver(h, albSrc,albDir, albIrr);
        end
        
end

function [varargout] = fastQuiver(h, xyz,uvw, tipsize)
        
        %vector consists of two lines, the butt and length of tip =tipsize*longest
        start = xyz';
        mid = (xyz .+ ((1.-tipsize).*uvw))';
        stop = (xyz + uvw)';
        blank = NaN( 1, size( xyz, 1));
        
        ubutt = [start(1,:); mid(1,:); blank](:);
        vbutt = [start(2,:); mid(2,:); blank](:);
        wbutt = [start(3,:); mid(3,:); blank](:);
        utipp = [mid(1,:); stop(1,:); blank](:);
        vtipp = [mid(2,:); stop(2,:); blank](:);
        wtipp = [mid(3,:); stop(3,:); blank](:);
        
        
        if !isequal(numel(h),2) || ...
           (!ishandle( h(1) ) || !ishandle( h(2)) ) ||...
           !isequal( [get( h, 'type' ){:}], "lineline" )
                
                h(2) = hbutt = line(ubutt,vbutt,wbutt);
                h(1) = htipp = line(utipp,vtipp,wtipp);
                
        else
                
                hbutt = h(2);
                set( hbutt,...
                        'xdata', ubutt,...
                        'ydata', vbutt,...
                        'zdata', wbutt ...
                );
                htipp = h(1);
                set( htipp,...
                        'xdata', utipp,...
                        'ydata', vtipp,...
                        'zdata', wtipp ...
                );
        end
        
        switch(nargout)
                case 0
                        return
                case 1
                        varargout{1} = h;
                case 2 
                        varargout{1} = htipp;
                        varargout{2} = hbutt;
        end
        
        nargoutchk(0,2);
    
end





%
%        SquareVertices=[0,1,1,0,0,1,1,0;...
%                        0,0,1,1,0,0,1,1;...
%                        0,0,0,0,1,1,1,1]';
%        SquareFaces=[1,2,6,5;...
%                     2,3,7,6;...
%                     3,4,8,7;...
%                     4,1,5,8;...
%                     1,2,3,4;...
%                     5,6,7,8];
%




