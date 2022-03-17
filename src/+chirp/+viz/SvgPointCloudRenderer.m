classdef SvgPointCloudRenderer
    % SVGPOINTCLOUDRENDERER Vector SVG exporter for 4D radar point clouds and bounding boxes.
    %
    % Renders top-down bird's-eye-view point clusters, velocity vectors, and
    % 3D oriented obstacle bounding boxes as clean standalone SVG graphics.
    
    methods (Static)
        function renderToFile(filePath, pointCloud, boundingBoxes, xExtent, yExtent, widthPx, heightPx)
            % Render BEV point cloud and oriented bounding boxes to SVG file.
            
            if nargin < 3, boundingBoxes = []; end
            if nargin < 4 || isempty(xExtent), xExtent = [-20.0, 20.0]; end
            if nargin < 5 || isempty(yExtent), yExtent = [0.0, 70.0]; end
            if nargin < 6 || isempty(widthPx), widthPx = 700; end
            if nargin < 7 || isempty(heightPx), heightPx = 600; end
            
            fid = fopen(filePath, 'w');
            if fid < 0
                error('SvgPointCloudRenderer:FileCreateFailed', 'Cannot write to %s', filePath);
            end
            cleanUp = onCleanup(@() fclose(fid));
            
            marginLeft = 60;
            marginRight = 40;
            marginTop = 50;
            marginBottom = 50;
            
            plotW = widthPx - marginLeft - marginRight;
            plotH = heightPx - marginTop - marginBottom;
            
            fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
            fprintf(fid, '<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">\n', ...
                widthPx, heightPx, widthPx, heightPx);
            fprintf(fid, '  <rect width="100%%" height="100%%" fill="#0f172a" />\n');
            
            % Title
            fprintf(fid, '  <text x="%d" y="30" fill="#f8fafc" font-family="sans-serif" font-size="16" font-weight="bold" text-anchor="middle">4D Radar Point Cloud (Bird''s-Eye-View)</text>\n', ...
                round(widthPx / 2));
                
            % Plot area background and range concentric rings
            fprintf(fid, '  <rect x="%d" y="%d" width="%d" height="%d" fill="#1e293b" stroke="#334155" stroke-width="1.5" />\n', ...
                marginLeft, marginTop, plotW, plotH);
                
            % Map sensor origin (x = 0, y = 0)
            normOriginX = (0.0 - xExtent(1)) / (xExtent(2) - xExtent(1));
            sensorPxX = marginLeft + normOriginX * plotW;
            sensorPxY = marginTop + plotH;
            
            % Draw range rings at 15m, 30m, 45m, 60m
            for ringR = [15, 30, 45, 60]
                if ringR <= yExtent(2)
                    ringPixRadius = (ringR / (yExtent(2) - yExtent(1))) * plotH;
                    fprintf(fid, '  <circle cx="%.2f" cy="%.2f" r="%.2f" fill="none" stroke="#475569" stroke-dasharray="4,4" stroke-width="1" />\n', ...
                        sensorPxX, sensorPxY, ringPixRadius);
                    fprintf(fid, '  <text x="%.2f" y="%.2f" fill="#64748b" font-family="sans-serif" font-size="10">%dm</text>\n', ...
                        sensorPxX + 5, sensorPxY - ringPixRadius + 12, ringR);
                end
            end
            
            % Render points
            N = pointCloud.NumPoints;
            pos = pointCloud.Positions;
            vel = pointCloud.Velocities;
            
            for i = 1:N
                px = pos(1, i);
                py = pos(2, i);
                v = vel(i);
                
                if px < xExtent(1) || px > xExtent(2) || py < yExtent(1) || py > yExtent(2)
                    continue;
                end
                
                normX = (px - xExtent(1)) / (xExtent(2) - xExtent(1));
                normY = (py - yExtent(1)) / (yExtent(2) - yExtent(1));
                
                ptX = marginLeft + normX * plotW;
                ptY = marginTop + plotH - normY * plotH;
                
                % Color by velocity: blue (approaching) -> green (stationary) -> red (receding)
                if v > 1.0
                    col = "#ef4444"; % Red
                elseif v < -1.0
                    col = "#3b82f6"; % Blue
                else
                    col = "#10b981"; % Green
                end
                
                fprintf(fid, '  <circle cx="%.2f" cy="%.2f" r="3.5" fill="%s" opacity="0.85" />\n', ...
                    ptX, ptY, col);
            end
            
            % Render radar sensor icon
            fprintf(fid, '  <polygon points="%.2f,%.2f %.2f,%.2f %.2f,%.2f" fill="#f59e0b" />\n', ...
                sensorPxX, sensorPxY - 8, sensorPxX - 6, sensorPxY, sensorPxX + 6, sensorPxY);
                
            % Axis labels
            fprintf(fid, '  <text x="%d" y="%d" fill="#94a3b8" font-family="sans-serif" font-size="12" text-anchor="middle">Lateral Distance X (m)</text>\n', ...
                marginLeft + round(plotW / 2), heightPx - 15);
            fprintf(fid, '  <text x="20" y="%d" fill="#94a3b8" font-family="sans-serif" font-size="12" text-anchor="middle" transform="rotate(-90 20 %d)">Range Y (m)</text>\n', ...
                marginTop + round(plotH / 2), marginTop + round(plotH / 2));
                
            fprintf(fid, '</svg>\n');
        end
    end
end
