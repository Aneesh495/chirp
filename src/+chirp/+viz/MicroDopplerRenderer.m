classdef MicroDopplerRenderer
    % MICRODOPPLERRENDERER Visualizer for micro-Doppler time-frequency spectrograms.
    %
    % Generates SVG vector spectrograms with overlaid kinematic cadence envelopes
    % and biomechanical stride indicators.
    
    methods (Static)
        function renderToFile(filePath, spectrogramDb, timeAxis, freqAxis, feats, widthPx, heightPx)
            % Render micro-Doppler spectrogram and tracking envelopes to SVG file.
            
            if nargin < 5, feats = []; end
            if nargin < 6 || isempty(widthPx), widthPx = 800; end
            if nargin < 7 || isempty(heightPx), heightPx = 450; end
            
            fid = fopen(filePath, 'w');
            if fid < 0
                error('MicroDopplerRenderer:FileCreateFailed', 'Cannot write to %s', filePath);
            end
            cleanUp = onCleanup(@() fclose(fid));
            
            [nF, nT] = size(spectrogramDb);
            maxVal = max(spectrogramDb(:));
            minVal = max(min(spectrogramDb(:)), maxVal - 50.0);
            
            marginLeft = 70;
            marginRight = 40;
            marginTop = 50;
            marginBottom = 50;
            
            plotW = widthPx - marginLeft - marginRight;
            plotH = heightPx - marginTop - marginBottom;
            
            fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
            fprintf(fid, '<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">\n', ...
                widthPx, heightPx, widthPx, heightPx);
            fprintf(fid, '  <rect width="100%%" height="100%%" fill="#090d16" />\n');
            
            fprintf(fid, '  <text x="%d" y="30" fill="#f1f5f9" font-family="sans-serif" font-size="16" font-weight="bold" text-anchor="middle">Micro-Doppler Time-Frequency Signature</text>\n', ...
                round(widthPx / 2));
                
            fprintf(fid, '  <rect x="%d" y="%d" width="%d" height="%d" fill="#0f172a" stroke="#334155" stroke-width="1.5" />\n', ...
                marginLeft, marginTop, plotW, plotH);
                
            % Resampled spectrogram rects
            stepF = max(1, floor(nF / 50));
            stepT = max(1, floor(nT / 80));
            cellW = plotW / ceil(nT / stepT);
            cellH = plotH / ceil(nF / stepF);
            
            subF = 1:stepF:nF;
            subT = 1:stepT:nT;
            
            for iF = 1:length(subF)
                fIdx = subF(iF);
                yPos = marginTop + plotH - iF * cellH;
                for iT = 1:length(subT)
                    tIdx = subT(iT);
                    xPos = marginLeft + (iT - 1) * cellW;
                    
                    val = spectrogramDb(fIdx, tIdx);
                    normVal = (val - minVal) / max(maxVal - minVal, eps);
                    normVal = max(0.0, min(1.0, normVal));
                    
                    if normVal > 0.05
                        % Inferno colormap approximation: dark purple -> magenta -> yellow
                        rCol = round(255 * min(1.0, 1.5 * normVal));
                        gCol = round(255 * max(0.0, 1.2 * normVal - 0.2));
                        bCol = round(255 * max(0.0, 0.8 * (1.0 - normVal)));
                        fprintf(fid, '  <rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" fill="rgb(%d,%d,%d)" />\n', ...
                            xPos, yPos, cellW + 0.5, cellH + 0.5, rCol, gCol, bCol);
                    end
                end
            end
            
            % Overlay envelope lines if features provided
            if ~isempty(feats) && isfield(feats, 'upperEnvelope')
                upEnv = feats.upperEnvelope;
                lowEnv = feats.lowerEnvelope;
                
                upPts = '';
                lowPts = '';
                for t = 1:length(upEnv)
                    xPix = marginLeft + (t / length(upEnv)) * plotW;
                    
                    normUp = (upEnv(t) - freqAxis(1)) / (freqAxis(end) - freqAxis(1));
                    normLow = (lowEnv(t) - freqAxis(1)) / (freqAxis(end) - freqAxis(1));
                    
                    yUp = marginTop + plotH - normUp * plotH;
                    yLow = marginTop + plotH - normLow * plotH;
                    
                    upPts = [upPts, sprintf('%.1f,%.1f ', xPix, yUp)];
                    lowPts = [lowPts, sprintf('%.1f,%.1f ', xPix, yLow)];
                end
                
                fprintf(fid, '  <polyline points="%s" fill="none" stroke="#38bdf8" stroke-width="2" />\n', strtrim(upPts));
                fprintf(fid, '  <polyline points="%s" fill="none" stroke="#38bdf8" stroke-width="2" />\n', strtrim(lowPts));
            end
            
            % Axis labels
            fprintf(fid, '  <text x="%d" y="%d" fill="#94a3b8" font-family="sans-serif" font-size="12" text-anchor="middle">Time (s)</text>\n', ...
                marginLeft + round(plotW / 2), heightPx - 15);
            fprintf(fid, '  <text x="20" y="%d" fill="#94a3b8" font-family="sans-serif" font-size="12" text-anchor="middle" transform="rotate(-90 20 %d)">Doppler Frequency (Hz)</text>\n', ...
                marginTop + round(plotH / 2), marginTop + round(plotH / 2));
                
            fprintf(fid, '</svg>\n');
        end
    end
end
