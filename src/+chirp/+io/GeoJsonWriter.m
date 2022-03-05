classdef GeoJsonWriter
    % GEOJSONWRITER Standard GeoJSON FeatureCollection exporter for radar perception data.
    %
    % Serializes 4D radar point clouds, 3D oriented bounding boxes, and tracked
    % trajectories into standard RFC 7946 GeoJSON format.
    
    methods (Static)
        function writePointCloud(filePath, pointCloud, originLatLon)
            % Export point cloud as GeoJSON FeatureCollection.
            %
            % Inputs:
            %   filePath     - Target output file path (.geojson)
            %   pointCloud   - chirp.pointcloud.PointCloud4D instance
            %   originLatLon - (Optional) [lat0, lon0] sensor origin coordinates
            
            if nargin < 3 || isempty(originLatLon)
                originLatLon = [37.7749, -122.4194]; % Reference coordinate
            end
            
            fid = fopen(filePath, 'w');
            if fid < 0
                error('GeoJsonWriter:FileCreateFailed', 'Cannot write to %s', filePath);
            end
            cleanUp = onCleanup(@() fclose(fid));
            
            lat0 = originLatLon(1);
            lon0 = originLatLon(2);
            
            % Conversion: meters to degrees lat/lon (equatorial approx)
            mPerDegLat = 111132.954;
            mPerDegLon = 111132.954 * cos(lat0 * pi / 180.0);
            
            N = pointCloud.NumPoints;
            pos = pointCloud.Positions;
            vel = pointCloud.Velocities;
            snr = pointCloud.SnrDb;
            
            fprintf(fid, '{\n  "type": "FeatureCollection",\n  "features": [\n');
            
            for i = 1:N
                lon = lon0 + (pos(1, i) / mPerDegLon);
                lat = lat0 + (pos(2, i) / mPerDegLat);
                alt = pos(3, i);
                
                fprintf(fid, '    {\n');
                fprintf(fid, '      "type": "Feature",\n');
                fprintf(fid, '      "geometry": {\n');
                fprintf(fid, '        "type": "Point",\n');
                fprintf(fid, '        "coordinates": [%.7f, %.7f, %.2f]\n', lon, lat, alt);
                fprintf(fid, '      },\n');
                fprintf(fid, '      "properties": {\n');
                fprintf(fid, '        "velocity": %.2f,\n', vel(i));
                fprintf(fid, '        "snr": %.2f\n', snr(i));
                fprintf(fid, '      }\n');
                if i < N
                    fprintf(fid, '    },\n');
                else
                    fprintf(fid, '    }\n');
                end
            end
            
            fprintf(fid, '  ]\n}\n');
        end
    end
end
