classdef RangeProcessor
    % RANGEPROCESSOR 1D fast-time Range-FFT engine with windowing and zero-padding.
    %
    % Transforms raw fast-time ADC beat samples into range profiles. Supports
    % DC offset removal, configurable windowing, zero-padding, and physical
    % range axis calibration.
    
    properties
        WindowType = 'blackman-harris'; % Fast-time window type
        Nfft = [];                      % FFT length (empty = next power of 2)
        RemoveDc = true;                % Subtract mean across fast-time
    end
    
    methods
        function obj = RangeProcessor(varargin)
            % Construct RangeProcessor with optional name-value configuration.
            if nargin > 0
                for k = 1:2:nargin
                    propName = varargin{k};
                    propVal = varargin{k+1};
                    if isprop(obj, propName)
                        obj.(propName) = propVal;
                    else
                        error('RangeProcessor:UnknownProperty', 'Property %s not found.', propName);
                    end
                end
            end
        end
        
        function [rangeProfiles, rangeAxis] = process(obj, radarCube, chirpObj)
            % Execute 1D Range-FFT across dimension 1 (fast-time) of the radar cube.
            %
            % Inputs:
            %   radarCube     - [N_adc x N_chirps x N_virtual] complex raw data
            %   chirpObj      - chirp.waveform.FmcwChirp configuration instance
            %
            % Returns:
            %   rangeProfiles - [N_range x N_chirps x N_virtual] complex spectra
            %   rangeAxis     - 1 x N_range physical range coordinates in meters
            
            [nAdc, nChirps, nVirt] = size(radarCube);
            
            % Determine FFT size
            if isempty(obj.Nfft)
                nFft = 2^chirp.utils.FastMath.nextpow2_fast(nAdc);
            else
                nFft = obj.Nfft;
            end
            
            % Generate window
            w = chirp.dsp.WindowFunction.create(obj.WindowType, nAdc);
            % Window normalization factor for amplitude preservation
            wNorm = w / sum(w);
            
            % Preallocate output
            rangeProfiles = zeros(nFft, nChirps, nVirt);
            
            % Process across virtual channels and chirps
            for iVirt = 1:nVirt
                sliceData = squeeze(radarCube(:, :, iVirt));
                
                % Remove DC bias if configured
                if obj.RemoveDc
                    sliceData = sliceData - mean(sliceData, 1);
                end
                
                % Apply fast-time window (broadcast across chirps)
                windowedData = sliceData .* wNorm;
                
                % Compute FFT along column dimension (fast-time)
                rangeProfiles(:, :, iVirt) = fft(windowedData, nFft, 1);
            end
            
            % Compute physical range axis
            % Beat frequency fb = k * (Fs / nFft)
            % Range R = (c * fb) / (2 * S)
            c = chirp.utils.PhysicalConstants.C;
            deltaF = chirpObj.SampleRate / nFft;
            kIndices = 0:(nFft - 1);
            beatFreqAxis = kIndices * deltaF;
            rangeAxis = (c * beatFreqAxis) / (2.0 * chirpObj.Slope);
        end
    end
end
