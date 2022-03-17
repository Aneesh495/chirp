classdef PolarPlotRenderer
    % POLARPLOTRENDERER Range-Azimuth Plan Position Indicator (PPI) radar scope renderer.
    %
    % Converts Range-Azimuth beamformed spatial power into polar sector displays
    % spanning -60 deg to +60 deg azimuth and 0 to R_max meters.
    
    methods (Static)
        function renderTerminal(rangeAzimuthDb, rangeAxis, azAxisDeg, outCols, outRows)
            % Render polar PPI radar scope in terminal.
            
            if nargin < 4 || isempty(outCols), outCols = 51; end
            if nargin < 5 || isempty(outRows), outRows = 21; end
            
            maxVal = max(rangeAzimuthDb(:));
            minVal = max(min(rangeAzimuthDb(:)), maxVal - 40.0);
            
            charGrid = repmat(' ', outRows, outCols);
            
            [nR, nAz] = size(rangeAzimuthDb);
            rMax = rangeAxis(end);
            
            % Center of radar origin at bottom center
            originCol = round(outCols / 2);
            originRow = outRows;
            
            for r = 1:outRows
                for c = 1:outCols
                    dx = (c - originCol);
                    dy = (originRow - r);
                    
                    if dy < 0
                        continue;
                    end
                    
                    distPix = sqrt(dx^2 + dy^2);
                    normDist = distPix / (outRows - 1);
                    
                    if normDist > 1.0
                        continue;
                    end
                    
                    angleRad = atan2(dx, dy);
                    angleDeg = angleRad * (180.0 / pi);
                    
                    if abs(angleDeg) > 60.0
                        continue;
                    end
                    
                    % Map to matrix coordinates
                    rIdx = max(1, min(nR, round(normDist * (nR - 1)) + 1));
                    azNorm = (angleDeg - azAxisDeg(1)) / (azAxisDeg(end) - azAxisDeg(1));
                    azIdx = max(1, min(nAz, round(azNorm * (nAz - 1)) + 1));
                    
                    val = rangeAzimuthDb(rIdx, azIdx);
                    normVal = (val - minVal) / max(maxVal - minVal, eps);
                    
                    if normVal > 0.8
                        charGrid(r, c) = '@';
                    elseif normVal > 0.6
                        charGrid(r, c) = '*';
                    elseif normVal > 0.4
                        charGrid(r, c) = ':';
                    elseif normVal > 0.2
                        charGrid(r, c) = '.';
                    end
                end
            end
            
            fprintf('\n=== Range-Azimuth Polar Scope (PPI) ===\n');
            fprintf('Azimuth: [%.1f .. %.1f] deg | Max Range: %.1f m\n', azAxisDeg(1), azAxisDeg(end), rMax);
            for r = 1:outRows
                fprintf('|%s|\n', charGrid(r, :));
            end
        end
    end
end
