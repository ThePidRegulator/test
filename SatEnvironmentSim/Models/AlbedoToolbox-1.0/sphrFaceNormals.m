


function faceNormals= sphrFaceNormals(map)
        CONST.EMR = 6371.01e3;

        [sy, sx] = size(map);

        faceNormals = zeros(sy, sx,3);

        [sphx, sphy, sphz] = spherMesh(sy+1,sx+1);
        faceNormals(:,:,1) = ( sphx(1:end-1, 1:end-1) + sphx(2:end, 1:end-1) +...
                              sphx(2:end, 2:end) + sphx(1:end-1, 2:end) )./4;

        faceNormals(:,:,2) = ( sphy(1:end-1, 1:end-1) + sphy(2:end, 1:end-1) +...
                              sphy(2:end, 2:end) + sphy(1:end-1, 2:end) )./4;

        faceNormals(:,:,3) = ( sphz(1:end-1, 1:end-1) + sphz(2:end, 1:end-1) +...
                              sphz(2:end, 2:end) + sphz(1:end-1, 2:end) )./4;

return
