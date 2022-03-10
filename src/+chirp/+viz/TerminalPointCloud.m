classdef TerminalPointCloud
    % TERMINALPOINTCLOUD Bird's-Eye-View (BEV) ANSI terminal radar point cloud renderer.
    %
    % Displays 2D top-down lateral vs longitudinal obstacle locations in the terminal,
    % coloring points by radial velocity.
    
    methods (Static)
        function render(pointCloud, xExtent, yExtent, outWidth, outHeight)
            % Print BEV scatter plot to terminal.
            %
            % Inputs:
            %   pointCloud - chirp.pointcloud.PointCloud4D instance
            %   xExtent    - [xMin, xMax] lateral display window (default: [-15, 15] m)
            %   yExtent    - [yMin, yMax] longitudinal display window (default: [0, 60] m)
            %   outWidth   - Character columns (default: 61)
            %   outHeight  - Character rows (default: 21)
            
            if nargin < 2 || isempty(xExtent), xExtent = [-15.0, 15.0]; end
            if nargin < 3 || isempty(yExtent), yExtent = [0.0, 60.0]; end
            if nargin < 4 || isempty(outWidth), outWidth = 61; end
            if nargin < 5 || isempty(outHeight), outHeight = 21; end
            
            % Grid buffer initialization
            charGrid = repmat(' ', outHeight, outWidth);
            
            N = pointCloud.NumPoints;
            pos = pointCloud.Positions;
            vel = pointCloud.Velocities;
            
            fprintf('\n=== Bird''s-Eye-View (BEV) Radar Point Cloud ===\n');
            fprintf('X (Lateral): [%.1f .. %.1f] m | Y (Range): [%.1f .. %.1f] m | Points: %d\n', ...
                xExtent(1), xExtent(2), yExtent(1), yExtent(2), N);
            fprintf('+%s+\n', repmat('-', 1, outWidth));
            
            for i = 1:N
                px = pos(1, i);
                py = pos(2, i);
                v = vel(i);
                
                if px < xExtent(1) || px > xExtent(2) || py < yExtent(1) || py > yExtent(2)
                    continue;
                end
                
                % Map coordinates to character row and column
                normX = (px - xExtent(1)) / (xExtent(2) - xExtent(1));
                normY = (py - yExtent(1)) / (yExtent(2) - yExtent(1));
                
                col = round(normX * (outWidth - 1)) + 1;
                % Invert row so range increases upwards
                row = outHeight - round(normY * (outHeight - 1));
                
                col = max(1, min(outWidth, col));
                row = max(1, min(outHeight, row));
                
                if v > 1.0
                    charGrid(row, col) = '^'; % Receding
                elseif v < -1.0
                    charGrid(row, col) = 'v'; % Approaching
                else
                    charGrid(row, col) = 'o'; % Stationary
                end
            end
            
            % Place sensor origin 'S' at center bottom
            colSensor = round(0.5 * (outWidth - 1)) + 1;
            charGrid(outHeight, colSensor) = 'S';
            
            for r = 1:outHeight
                fprintf('|%s|\n', charGrid(r, :));
            end
            
            fprintf('+%s+\n', repmat('-', 1, outWidth));
            fprintf('Legend: (o) Stationary  (^) Receding  (v) Approaching  (S) Radar Sensor\n');
        end
    end
end
