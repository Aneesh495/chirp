classdef TerminalHeatmap
    % TERMINALHEATMAP ANSI terminal color heatmap renderer for Range-Doppler and SAR imagery.
    %
    % Renders 2D matrix power distributions directly into the command terminal
    % using ANSI 256-color blocks without requiring a graphical display server.
    
    methods (Static)
        function render(matrixDb, rowAxis, colAxis, titleStr, outWidth, outHeight)
            % Print an ANSI colored heatmap to the terminal standard output.
            %
            % Inputs:
            %   matrixDb  - [N_rows x N_cols] 2D matrix in decibels
            %   rowAxis   - 1 x N_rows coordinate vector (e.g. Range in m)
            %   colAxis   - 1 x N_cols coordinate vector (e.g. Velocity in m/s)
            %   titleStr  - Header title string
            %   outWidth  - Terminal character column width (default: 60)
            %   outHeight - Terminal row height (default: 20)
            
            if nargin < 4 || isempty(titleStr), titleStr = 'Range-Doppler Power Spectrum'; end
            if nargin < 5 || isempty(outWidth), outWidth = 60; end
            if nargin < 6 || isempty(outHeight), outHeight = 20; end
            
            [nR, nC] = size(matrixDb);
            
            % Downsample / resample matrix to terminal dimensions
            rIndices = round(linspace(1, nR, outHeight));
            cIndices = round(linspace(1, nC, outWidth));
            sampled = matrixDb(rIndices, cIndices);
            
            maxVal = max(sampled(:));
            minVal = max(min(sampled(:)), maxVal - 50.0); % 50 dB dynamic range
            
            fprintf('\n=== %s ===\n', titleStr);
            fprintf('Range: [%.1f .. %.1f] m | Doppler: [%.1f .. %.1f] m/s | Dynamic: [%.1f .. %.1f] dB\n', ...
                rowAxis(1), rowAxis(end), colAxis(1), colAxis(end), minVal, maxVal);
            fprintf('+%s+\n', repmat('-', 1, outWidth));
            
            % ANSI 256-color spectrum: 16 (black) -> 17..21 (blue) -> 28..34 (green) -> 220..226 (yellow) -> 196 (red)
            for r = 1:outHeight
                rowStr = '';
                for c = 1:outWidth
                    val = sampled(r, c);
                    normVal = (val - minVal) / max(maxVal - minVal, eps);
                    normVal = max(0.0, min(1.0, normVal));
                    
                    % Map normalized [0, 1] to ANSI 256 color code
                    if normVal < 0.15
                        colorCode = 233; % Very dark gray
                        ch = ' ';
                    elseif normVal < 0.35
                        colorCode = 24;  % Deep blue
                        ch = '.';
                    elseif normVal < 0.55
                        colorCode = 34;  % Green
                        ch = ':';
                    elseif normVal < 0.75
                        colorCode = 220; % Yellow
                        ch = '*';
                    elseif normVal < 0.90
                        colorCode = 208; % Orange
                        ch = '#';
                    else
                        colorCode = 196; % Bright red
                        ch = '@';
                    end
                    
                    % ANSI escape: \033[38;5;{colorCode}m{ch}\033[0m
                    rowStr = [rowStr, sprintf('\033[38;5;%dm%c', colorCode, ch)];
                end
                rowStr = [rowStr, sprintf('\033[0m')];
                fprintf('|%s|\n', rowStr);
            end
            
            fprintf('+%s+\n', repmat('-', 1, outWidth));
        end
    end
end
