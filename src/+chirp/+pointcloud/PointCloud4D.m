classdef PointCloud4D
    % POINTCLOUD4D Structured container and transformer for 4D kinematic point clouds.
    %
    % Manages high-density point buffers with Cartesian positions [3 x N],
    % radial velocities [1 x N], SNR values [1 x N], and RCS estimates [1 x N].
    
    properties
        Positions = zeros(3, 0);      % 3 x N matrix of Cartesian points [x; y; z] (m)
        Velocities = zeros(1, 0);     % 1 x N vector of radial velocities (m/s)
        SnrDb = zeros(1, 0);          % 1 x N vector of Signal-to-Noise Ratios (dB)
        Rcs = zeros(1, 0);            % 1 x N vector of Radar Cross Sections (m^2)
        Timestamp = 0.0;              % Acquisition frame epoch timestamp (s)
    end
    
    properties (Dependent)
        NumPoints                     % Number of active points in cloud N
    end
    
    methods
        function obj = PointCloud4D(pos, vel, snr, rcs, ts)
            % Construct PointCloud4D from component arrays.
            if nargin >= 1 && ~isempty(pos)
                if size(pos, 1) ~= 3
                    error('PointCloud4D:InvalidDimension', 'Positions must have 3 rows [x; y; z].');
                end
                obj.Positions = pos;
                n = size(pos, 2);
                
                if nargin >= 2 && ~isempty(vel), obj.Velocities = vel(:).'; else, obj.Velocities = zeros(1, n); end
                if nargin >= 3 && ~isempty(snr), obj.SnrDb = snr(:).'; else, obj.SnrDb = zeros(1, n); end
                if nargin >= 4 && ~isempty(rcs), obj.Rcs = rcs(:).'; else, obj.Rcs = ones(1, n); end
                if nargin >= 5 && ~isempty(ts), obj.Timestamp = ts; end
            end
        end
        
        function n = get.NumPoints(obj)
            n = size(obj.Positions, 2);
        end
        
        function cropped = cropBox(obj, xMin, xMax, yMin, yMax, zMin, zMax)
            % Filter points lying strictly within 3D bounding box boundaries.
            if obj.NumPoints == 0
                cropped = obj;
                return;
            end
            
            mask = (obj.Positions(1, :) >= xMin) & (obj.Positions(1, :) <= xMax) & ...
                   (obj.Positions(2, :) >= yMin) & (obj.Positions(2, :) <= yMax) & ...
                   (obj.Positions(3, :) >= zMin) & (obj.Positions(3, :) <= zMax);
                   
            cropped = chirp.pointcloud.PointCloud4D(...
                obj.Positions(:, mask), ...
                obj.Velocities(mask), ...
                obj.SnrDb(mask), ...
                obj.Rcs(mask), ...
                obj.Timestamp ...
            );
        end
        
        function filtered = filterVelocity(obj, vMin, vMax)
            % Filter points within radial velocity interval [vMin, vMax].
            if obj.NumPoints == 0
                filtered = obj;
                return;
            end
            
            mask = (obj.Velocities >= vMin) & (obj.Velocities <= vMax);
            filtered = chirp.pointcloud.PointCloud4D(...
                obj.Positions(:, mask), ...
                obj.Velocities(mask), ...
                obj.SnrDb(mask), ...
                obj.Rcs(mask), ...
                obj.Timestamp ...
            );
        end
        
        function transformed = transform(obj, rotMatrix, translationVec)
            % Apply 3D rigid Euclidean transformation: p_new = R * p + t.
            if obj.NumPoints == 0
                transformed = obj;
                return;
            end
            
            t = translationVec(:);
            newPos = rotMatrix * obj.Positions + t;
            
            transformed = chirp.pointcloud.PointCloud4D(...
                newPos, ...
                obj.Velocities, ...
                obj.SnrDb, ...
                obj.Rcs, ...
                obj.Timestamp ...
            );
        end
        
        function merged = append(obj, other)
            % Concatenate two point clouds into a unified container.
            merged = chirp.pointcloud.PointCloud4D(...
                [obj.Positions, other.Positions], ...
                [obj.Velocities, other.Velocities], ...
                [obj.SnrDb, other.SnrDb], ...
                [obj.Rcs, other.Rcs], ...
                obj.Timestamp ...
            );
        end
    end
end
