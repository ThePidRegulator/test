% function Vectors = CreateAlbedoVectors(rsat_e,Albedo_e)
%
% This function is based on 
% function P = ss_proj(re,n,a);
% re = satelite nadir vector in ECEF
% n = normal in ECEF that albedo contributions are projected onto 
% a = albedo matrix
% Made bye Dan D.V. Bhanderi
%
% Because Albedo vectors was more practical to use in the context of a sun
% sensor using several  photodiodes ss_proj was rewritten and renamed for 
% this purpose. 
% 
% Input:
% Albedo_e = is the Albedo matrix used to create the grid vector in ECEF
%            frame
% rsat_e = satelite position in ECEF frame
%
% Output:
% Vectors(1:3,:) = Albedo Vectors pointing from grid to satelite
% Vectors(4,:) = Reflected light from earth with sun intensity = CONST.AM0
% 
% Internal variables
% [grid(1), grid(2), grid(3)] = Grid vector giving grid position in ECEF, 
%                               points from center 
%
% Created: 07.04.2015, Henrik Rudi Haave
% 

function Vectors = CreateAlbedoVectors(rsat_e,Albedo_e)

CONST.EMR = 6371.01e3;
CONST.AM0 = 1366.9;

grid = zeros(3,1);

[sy, sx] = size(Albedo_e);

Vect = 0;
Vectors = zeros(4,300); 

for i=1:sy
  for j=1:sx
    if Albedo_e(i,j) > 0
      Vect = Vect +1;
      % Grid vector in ECEF
      [grid_theta, grid_phi] = idx2rad(i,j,sy,sx);
      [grid(1), grid(2), grid(3)] = sph2cart(grid_theta,pi/2-grid_phi,CONST.EMR);
      Vectors(4,Vect) = Albedo_e(i,j);
      Vectors(1:3,Vect) = (-rsat_e+grid)/norm(-rsat_e+grid);
      
    end
  end
end

return
