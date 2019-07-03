
% quiver3 is exrutiatingly slow and consists of 3 line elements. vectors are
% therefor made of Two lines, one short and one long in diferent colors or style to indicate direction

function IllustratetOrbitEnvironmentSim(steps,refl)
	
        
        EMR = 6371.010e3;
        AM0 = 1366.9;
        d2r = pi/180;
        AU2m = 149597870700;
                
        illuData = IllustrationData;
        % sy row
        % sx col
        %% Format and prep the data neded for illustration, timesteps are stored in {:}
        albedo = [steps(:).a_sat_ecef];
        mapSize = size( refl.data );
        
        % max number of steps, satelite and sun position
        illuData.S = size( steps, 2 );
        illuData.rsat = {steps(:).r_sat_ecef};
        illuData.rsun = {steps(:).r_sun_ecef};
        % Albedo directions, strengths and max reflected sunlight
        illuData.albDir = {albedo(:).vect};
        illuData.albIrr = {albedo(:).irr};
        illuData.maxIrr = max( ...
                [cellfun( @max, ...
                        {albedo(:).irr}, ...
                        'UniformOutput', false ...
                ){:}]
        );
        
        satsph = cart2sph( [illuData.rsat{:}]' );
        sunsph = cart2sph( [illuData.rsun{:}]' );

        % Convert phi to polar angle
        satsph(:,2) = pi/2 - satsph(:,2);
        sunsph(:,2) = pi/2 - sunsph(:,2);

        satFov = earthfov( satsph, refl );
        sunlit = earthfov( sunsph, refl );
        litFov = sunlit & satFov;
        unionIdx =find( litFov );
        
        
        % Cartesean Earth sources and indexes of wich sunlight is reflected from
        illuData.cartSrcMap = refl.gridPosMap;
        illuData.albSrcIdx = ...
        cellfun( @repmat, ...
                {albedo(:).unionMap}, ...
                {[1,1,3]}, ...
                'uniformOutput', false...
        );
       
      % Earth map of reflectivity masked by Sunlit area and satelite FOV
        sunlitRefl = mask( refl.data, unionMap, 0.2 );
        illuData.sunlitRefl = num2cell( sunlitRefl, [1, 2] );
        fovSunlitRefl = mask( sunlitRefl, sunMap, 0.5 );
        illuData.fovSunlitRefl = num2cell( fovSunlitRefl, [1, 2] );
        
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
        
%         % Illustration controll
%         sliderH = uicontrol (gcf,                    ...
%                'style', 'slider',                ...
%                'Units', 'normalized',            ...
%                'position', [0.5, 0.1, 0.8, 0.1], ...
%                'min', 1,                         ...
%                'max', IlluData.S,                        ...
%                'value', IlluData.s,                      ...
%                'callback', {@sliderStepping,illuData},...
%                'SliderStep', [1/(S-1) , 1/(S-1) ]...
%        );
        
        for s=2:25
                illuSteppTo(s, illuData)
                drawnow('expose');
                pause(0.1);
        end
        

        
        
        
%        clear "illuSatOrbit"
%        clear "illuSunIrr"
%        clear "illuAlbedoIrr"
%        clf
         %camlookat([rsunH,rsatH])
end



function sliderStepping(hslider ,~ , illuData)
        
        % Get step selected by slider
        s = ceil(get( hslider,"value" ))
        
        illustrationSteppTo(s, illuData)

end

function illuSteppTo(s, illuData)
        
        
        EMR = 6371.01e3;
        AM0 = 1366.9;
        d2r = pi/180;
        AU2m = 149597870700;
        
        rsat = illuData.rsat{s};
        
        scaling = norm(rsat) - EMR;
        
        illuSatOrbit(s, illuData.rsat);
     
        illuSunIrr(rsat, illuData.rsun{s}, scaling);
        
        albedoIrr = illuData.albIrr{s};
        maxIrr = illuData.maxIrr;
        albedoDirr = illuData.albDir{s}*scaling;
        mapIdx = illuData.albSrcIdx{s};
        albedoSrce = reshape( illuData.cartSrcMap(mapIdx), [], 3 );
        
        tic
        illuAlbedoIrr( albedoSrce, albedoDirr, albedoIrr ,maxIrr )
        toc
        
        set( illuData.EarthH,'Cdata', illuData.sunlitRefl{s}' );
       
end


function illuSatOrbit(s, rsat)
        
        stail = 50;
        
        persistent h;
        if isempty(h) || !ishandle(h)
                h = line( rsat{1}(1), rsat{1}(2), rsat{1}(3) );
                set( h,...
                        'linewidth', 2,...
                        'color', [0,1,0]...
                );
        else
                stail=s-stail;
                if( stail<1 )
                        stail=1;
                end
                
                rsat = reshape( [rsat{stail:s}], 3, [] );
                
                set( h,...
                        "xdata", rsat(1, :),...
                        "ydata", rsat(2, :),...
                        "zdata", rsat(3, :)...
                );
        end
end
        


function illuSunIrr(pointTo, rsun, scaling)
        AU2m = 149597870700;
        AM0 = 1366.9;
        
        sun = AM0 * 2000 * -rsun' / AU2m; % rsat points to the source, reverse it to make it point to what is iluminates
        sunxyz = pointTo' - sun;
        
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

function illuAlbedoIrr(albSrc, albDir, albIrr, maxIrr )
        
        albIrr = albIrr/max(albIrr);
        
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




