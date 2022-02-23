classdef Dca1000Reader
    % DCA1000READER Ingestion parser for Texas Instruments DCA1000 EVM raw LVDS captures.
    %
    % Unpacks 16-bit signed integer interleaved I/Q samples collected from
    % TI mmWave radar sensors (e.g. AWR1642, AWR1843, AWR2243) into 3D radar cubes.
    
    properties
        NumRx = 4;                    % Number of physical RX channels
        NumAdcSamples = 256;          % Fast-time ADC samples per chirp
        NumChirpsPerFrame = 64;       % Slow-time chirps per frame
        IsComplex = true;             % Complex 1x or 2x mode
    end
    
    methods
        function obj = Dca1000Reader(varargin)
            % Construct a Dca1000Reader instance.
            if nargin > 0
                for k = 1:2:nargin
                    propName = varargin{k};
                    propVal = varargin{k+1};
                    if isprop(obj, propName)
                        obj.(propName) = propVal;
                    else
                        error('Dca1000Reader:UnknownProperty', 'Property %s not found.', propName);
                    end
                end
            end
        end
        
        function radarCube = readBinaryFile(obj, filePath)
            % Read and reshape raw binary .bin file into 3D radar cube.
            %
            % Returns:
            %   radarCube - [NumAdcSamples x NumChirpsPerFrame x NumRx] complex double array
            
            fid = fopen(filePath, 'rb');
            if fid < 0
                error('Dca1000Reader:FileNotFound', 'Could not open binary capture file: %s', filePath);
            end
            
            cleanUp = onCleanup(@() fclose(fid));
            
            % Read 16-bit signed integers (int16)
            rawInts = fread(fid, Inf, 'int16');
            
            samplesPerChirp = obj.NumAdcSamples * obj.NumRx * 2; % 2 for I and Q
            totalChirps = floor(length(rawInts) / samplesPerChirp);
            
            if totalChirps < obj.NumChirpsPerFrame
                warning('Dca1000Reader:TruncatedFile', ...
                    'File contains %d chirps, fewer than configured frame size %d.', ...
                    totalChirps, obj.NumChirpsPerFrame);
                nChirps = totalChirps;
            else
                nChirps = obj.NumChirpsPerFrame;
            end
            
            % Truncate to exact frame size
            neededInts = obj.NumAdcSamples * obj.NumRx * 2 * nChirps;
            frameInts = rawInts(1:neededInts);
            
            % Reshape: 2 (I/Q) x NumRx x NumAdcSamples x nChirps
            % In standard DCA1000 EVM format:
            % samples are arranged as: [I_rx1, Q_rx1, I_rx2, Q_rx2, ...]
            reshaped = reshape(frameInts, [2, obj.NumRx, obj.NumAdcSamples, nChirps]);
            
            % Form complex double: I + 1j * Q
            complexData = squeeze(reshaped(1, :, :, :) + 1j * reshaped(2, :, :, :));
            
            % Permute dimensions to standard [NumAdcSamples x nChirps x NumRx]
            % Current order: [NumRx, NumAdcSamples, nChirps]
            radarCube = permute(complexData, [2, 3, 1]);
        end
        
        function writeSyntheticBinary(obj, radarCube, filePath)
            % Write a 3D radar cube into standard DCA1000 EVM binary format.
            
            [nAdc, nChirps, nRx] = size(radarCube);
            
            % Permute to [NumRx, NumAdcSamples, nChirps]
            permuted = permute(radarCube, [3, 1, 2]);
            
            % Interleave I and Q
            interleaved = zeros(2, nRx, nAdc, nChirps);
            interleaved(1, :, :, :) = real(permuted);
            interleaved(2, :, :, :) = imag(permuted);
            
            % Scale to int16 range
            maxVal = max(abs(radarCube(:)));
            if maxVal > 0
                scale = 32767.0 / maxVal;
            else
                scale = 1.0;
            end
            
            intData = int16(round(interleaved(:) * scale));
            
            fid = fopen(filePath, 'wb');
            if fid < 0
                error('Dca1000Reader:FileCreateFailed', 'Cannot write to file: %s', filePath);
            end
            cleanUp = onCleanup(@() fclose(fid));
            
            fwrite(fid, intData, 'int16');
        end
    end
end
