function [x,y,z] = spherMesh( m, n )
        
        dn = pi*2/n;
        dm = pi/m;
        
        theta = linspace (-pi , pi, n);
        phi = linspace (-pi/2 , pi/2, m);
        [theta,phi] = meshgrid (theta, phi);
        
        x = cos (phi) .* cos (theta);
        y = cos (phi) .* sin (theta);
        z = sin (phi);
	
end