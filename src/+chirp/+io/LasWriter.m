classdef LasWriter
    % LASWRITER ASPRS LAS 1.2 binary point cloud format exporter.
    %
    % Writes standard geospatial binary LAS point files containing scaled XYZ coordinates,
    % integer intensity, scan angle, and classification codes.
    
    methods (Static)
        function write(filePath, pointCloud)
            % Write PointCloud4D to standard binary LAS 1.2 file format.
            
            fid = fopen(filePath, 'wb');
            if fid < 0
                error('LasWriter:FileCreateFailed', 'Cannot write to %s', filePath);
            end
            cleanUp = onCleanup(@() fclose(fid));
            
            N = pointCloud.NumPoints;
            
            % LAS 1.2 Header (227 bytes)
            fwrite(fid, 'LASF', 'char');          % File Signature (4 bytes)
            fwrite(fid, 0, 'uint16');             % File Source ID
            fwrite(fid, 0, 'uint16');             % Global Encoding
            fwrite(fid, zeros(4, 1), 'uint32');   % Project ID GUID
            fwrite(fid, 1, 'uint8');              % Version Major (1)
            fwrite(fid, 2, 'uint8');              % Version Minor (2)
            fwrite(fid, repmat(' ', 1, 32), 'char'); % System Identifier
            fwrite(fid, repmat(' ', 1, 32), 'char'); % Generating Software
            fwrite(fid, 0, 'uint16');             % Flight Date DOY
            fwrite(fid, 2022, 'uint16');          % Year
            fwrite(fid, 227, 'uint16');           % Header Size
            fwrite(fid, 227, 'uint32');           % Offset to Point Data
            fwrite(fid, 0, 'uint32');             % Number of Variable Length Records
            fwrite(fid, 1, 'uint8');              % Point Data Format (Format 1: XYZ, Intensity, GPS time)
            fwrite(fid, 28, 'uint16');            % Point Data Record Length (28 bytes for Format 1)
            fwrite(fid, N, 'uint32');             % Legacy Number of Point Records
            fwrite(fid, zeros(5, 1), 'uint32');   % Number of points by return
            
            % Scale factors and offsets (0.001 m resolution)
            scaleXYZ = [0.001; 0.001; 0.001];
            fwrite(fid, scaleXYZ, 'double');      % X, Y, Z scale factors
            fwrite(fid, zeros(3, 1), 'double');   % X, Y, Z offsets
            
            if N > 0
                minPos = min(pointCloud.Positions, [], 2);
                maxPos = max(pointCloud.Positions, [], 2);
            else
                minPos = zeros(3, 1);
                maxPos = zeros(3, 1);
            end
            
            fwrite(fid, maxPos(1), 'double');     % Max X
            fwrite(fid, minPos(1), 'double');     % Min X
            fwrite(fid, maxPos(2), 'double');     % Max Y
            fwrite(fid, minPos(2), 'double');     % Min Y
            fwrite(fid, maxPos(3), 'double');     % Max Z
            fwrite(fid, minPos(3), 'double');     % Min Z
            
            % Write Point Records (Format 1: 28 bytes per point)
            pos = pointCloud.Positions;
            snr = pointCloud.SnrDb;
            ts = pointCloud.Timestamp;
            
            for i = 1:N
                % Scaled integer coordinates
                xi = int32(round(pos(1, i) / 0.001));
                yi = int32(round(pos(2, i) / 0.001));
                zi = int32(round(pos(3, i) / 0.001));
                
                % Intensity (scale SNR 0..40 dB to uint16 0..65535)
                intensity = uint16(min(65535, max(0, round(snr(i) * 1600.0))));
                
                fwrite(fid, xi, 'int32');
                fwrite(fid, yi, 'int32');
                fwrite(fid, zi, 'int32');
                fwrite(fid, intensity, 'uint16');
                fwrite(fid, 1, 'uint8');          % Return number: 1
                fwrite(fid, 1, 'uint8');          % Classification: 1 (Unclassified)
                fwrite(fid, 0, 'int8');           % Scan angle
                fwrite(fid, 0, 'uint8');          % User data
                fwrite(fid, 0, 'uint16');         % Point source ID
                fwrite(fid, ts, 'double');        % GPS Time
            end
        end
    end
end
