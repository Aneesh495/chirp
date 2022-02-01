classdef StftAnalyzer
    % STFTANALYZER Short-Time Fourier Transform time-frequency spectrogram generator.
    %
    % Computes discrete STFT and calibrated micro-Doppler spectrograms from
    % complex slow-time radar pulse returns.
    
    properties
        WindowLength = 64;            % Window length in samples
        HopSize = 8;                  % Step size between successive windows (overlap = WindowLength - HopSize)
        Nfft = 128;                   % FFT evaluation length
        WindowType = 'blackman-harris'; % Window function type
    end
    
    methods
        function obj = StftAnalyzer(varargin)
            % Construct an StftAnalyzer with optional name-value configuration.
            if nargin > 0
                for k = 1:2:nargin
                    propName = varargin{k};
                    propVal = varargin{k+1};
                    if isprop(obj, propName)
                        obj.(propName) = propVal;
                    else
                        error('StftAnalyzer:UnknownProperty', 'Property %s not found.', propName);
                    end
                end
            end
        end
        
        function [spectrogramDb, timeAxis, freqAxis] = compute(obj, signal, fs)
            % Compute centered STFT spectrogram of a 1D complex signal.
            %
            % Inputs:
            %   signal         - 1D complex signal vector (e.g. slow-time phase returns)
            %   fs             - Sampling rate in Hz (e.g. PRF)
            %
            % Returns:
            %   spectrogramDb  - [N_freq x N_time] power spectrogram in decibels
            %   timeAxis       - 1 x N_time center time vector (s)
            %   freqAxis       - 1 x N_freq centered Doppler frequency vector (Hz)
            
            x = signal(:).';
            N = length(x);
            winLen = obj.WindowLength;
            hop = obj.HopSize;
            nFft = obj.Nfft;
            
            if winLen > N
                winLen = N;
                obj.WindowLength = winLen;
            end
            if nFft < winLen
                nFft = 2^chirp.utils.FastMath.nextpow2_fast(winLen);
                obj.Nfft = nFft;
            end
            
            % Generate window
            w = chirp.dsp.WindowFunction.create(obj.WindowType, winLen);
            w = (w / sum(w)).';
            
            % Compute number of frames
            numFrames = floor((N - winLen) / hop) + 1;
            stftMatrix = zeros(nFft, numFrames);
            timeIndices = zeros(1, numFrames);
            
            for k = 1:numFrames
                startIdx = (k - 1) * hop + 1;
                segment = x(startIdx:(startIdx + winLen - 1)) .* w;
                
                % Centered FFT
                fftCol = fftshift(fft(segment, nFft));
                stftMatrix(:, k) = fftCol(:);
                timeIndices(k) = (startIdx + (winLen - 1) / 2.0) / fs;
            end
            
            pwrLin = abs(stftMatrix).^2;
            spectrogramDb = 10.0 * log10(max(pwrLin, eps));
            
            timeAxis = timeIndices;
            
            if mod(nFft, 2) == 0
                kIndices = (-nFft/2):(nFft/2 - 1);
            else
                kIndices = (-(nFft - 1)/2):((nFft - 1)/2);
            end
            freqAxis = kIndices * (fs / nFft);
        end
    end
end
