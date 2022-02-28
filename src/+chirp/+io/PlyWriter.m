classdef PlyWriter
    % PLYWRITER Stanford Triangle Format (.ply) 3D point cloud exporter.
    %
    % Exports 4D radar point clouds with Cartesian coordinates (x, y, z),
    % Doppler velocity mapped to RGB colormaps, and SNR intensity attributes.
    
    methods (Static)
        function writeAscii(filePath, pointCloud)
            % Write point cloud to ASCII .ply file.
            
            fid = fopen(filePath, 'w');
            if fid < 0
                error('PlyWriter:FileCreateFailed', 'Cannot write to %s', filePath);
            end
            cleanUp = onCleanup(@() fclose(fid));
            
            N = pointCloud.NumPoints;
            
            % PLY Header
            fprintf(fid, 'ply\n');
            fprintf(fid, 'format ascii 1.0\n');
            fprintf(fid, 'comment Created by chirp 4D radar engine\n');
            fprintf(fid, 'element vertex %d\n', N);
            fprintf(fid, 'property float x\n');
            fprintf(fid, 'property float y\n');
            fprintf(fid, 'property float z\n');
            fprintf(fid, 'property float velocity\n');
            fprintf(fid, 'property float snr\n');
            fprintf(fid, 'property uchar red\n');
            fprintf(fid, 'property uchar green\n');
            fprintf(fid, 'property uchar blue\n');
            fprintf(fid, 'end_header\n');
            
            if N == 0
                return;
            end
            
            pos = pointCloud.Positions;
            vel = pointCloud.Velocities;
            snr = pointCloud.SnrDb;
            
            % Map velocities to Jet RGB colors
            vMin = -20.0; vMax = 20.0;
            normVel = (vel - vMin) / (vMax - vMin);
            normVel = max(0.0, min(1.0, normVel));
            
            for i = 1:N
                nv = normVel(i);
                % Jet color mapping: blue -> cyan -> green -> yellow -> red
                rVal = uint8(255 * max(0.0, min(1.0, 1.5 - abs(4.0 * nv - 3.0))));
                gVal = uint8(255 * max(0.0, min(1.0, 1.5 - abs(4.0 * nv - 2.0))));
                bVal = uint8(255 * max(0.0, min(1.0, 1.5 - abs(4.0 * nv - 1.0))));
                
                fprintf(fid, '%.4f %.4f %.4f %.3f %.2f %d %d %d\n', ...
                    pos(1, i), pos(2, i), pos(3, i), vel(i), snr(i), rVal, gVal, bVal);
            end
        end
    end
end
