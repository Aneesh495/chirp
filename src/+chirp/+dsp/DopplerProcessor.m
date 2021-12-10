classdef DopplerProcessor
    % DOPPLERPROCESSOR 2D slow-time Doppler-FFT processor and velocity mapper.
    %
    % Transforms Range-Chirp profiles into Range-Doppler matrices (RDM) across
    % slow-time pulse dimension. Applies windowing, zero-padding, fftshift,
    % and physical velocity axis mapping.
    
    properties
        WindowType = 'chebyshev';       % Slow-time window type
        ChebyshevAttenuationDb = 60.0; % Chebyshev sidelobe attenuation
        Nfft = [];                     % Doppler FFT size (empty = next power of 2)
    end
    
    methods
        function obj = DopplerProcessor(varargin)
            % Construct DopplerProcessor with optional name-value overrides.
            if nargin > 0
                for k = 1:2:nargin
                    propName = varargin{k};
                    propVal = varargin{k+1};
                    if isprop(obj, propName)
                        obj.(propName) = propVal;
                    else
                        error('DopplerProcessor:UnknownProperty', 'Property %s not found.', propName);
                    end
                end
            end
        end
        
        function [rdmCube, velocityAxis] = process(obj, rangeProfiles, chirpObj)
            % Execute slow-time Doppler-FFT along dimension 2 of rangeProfiles.
            %
            % Inputs:
            %   rangeProfiles - [N_range x N_chirps x N_virtual] complex data
            %   chirpObj      - chirp.waveform.FmcwChirp configuration instance
            %
            % Returns:
            %   rdmCube       - [N_range x N_doppler x N_virtual] centered RDM cube
            %   velocityAxis  - 1 x N_doppler physical radial velocity axis (m/s)
            
            [nRange, nChirps, nVirt] = size(rangeProfiles);
            
            if isempty(obj.Nfft)
                nFft = 2^chirp.utils.FastMath.nextpow2_fast(nChirps);
            else
                nFft = obj.Nfft;
            end
            
            % Generate slow-time window
            w = chirp.dsp.WindowFunction.create(obj.WindowType, nChirps, obj.ChebyshevAttenuationDb);
            wNorm = (w / sum(w)).';
            
            rdmCube = zeros(nRange, nFft, nVirt);
            
            for iVirt = 1:nVirt
                sliceData = squeeze(rangeProfiles(:, :, iVirt));
                
                % Window across slow-time chirps
                windowedData = sliceData .* wNorm;
                
                % FFT along dimension 2 (chirps), followed by fftshift to center 0 m/s
                fftData = fft(windowedData, nFft, 2);
                rdmCube(:, :, iVirt) = fftshift(fftData, 2);
            end
            
            % Compute physical velocity axis
            % PRF = 1 / T_chirp
            % delta_f_doppler = PRF / nFft
            % v = (lambda * f_doppler) / 2
            lambda = chirpObj.Wavelength;
            prf = 1.0 / chirpObj.TotalChirpPeriod;
            
            if mod(nFft, 2) == 0
                kIndices = (-nFft/2):(nFft/2 - 1);
            else
                kIndices = (-(nFft - 1)/2):((nFft - 1)/2);
            end
            
            dopplerFreqAxis = kIndices * (prf / nFft);
            % Doppler formula: fd = 2 * v / lambda => v = fd * lambda / 2
            velocityAxis = dopplerFreqAxis * (lambda / 2.0);
        end
    end
end
