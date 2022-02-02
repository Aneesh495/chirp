classdef SynchrosqueezingTransform
    % SYNCHROSQUEEZINGTRANSFORM Fourier Synchrosqueezed time-frequency reallocation.
    %
    % Sharpens diffuse STFT spectrogram ridges by reallocating energy along the
    % instantaneous frequency axis, concentrating kinematic Doppler modulations.
    
    properties
        WindowLength = 64;
        HopSize = 4;
        Nfft = 128;
    end
    
    methods
        function obj = SynchrosqueezingTransform(varargin)
            % Construct an instance with optional parameters.
            if nargin > 0
                for k = 1:2:nargin
                    propName = varargin{k};
                    propVal = varargin{k+1};
                    if isprop(obj, propName)
                        obj.(propName) = propVal;
                    else
                        error('SynchrosqueezingTransform:UnknownProperty', 'Property %s not found.', propName);
                    end
                end
            end
        end
        
        function [sstMatrix, timeAxis, freqAxis] = compute(obj, signal, fs)
            % Compute Fourier Synchrosqueezed Transform of 1D signal.
            
            x = signal(:).';
            N = length(x);
            winLen = obj.WindowLength;
            hop = obj.HopSize;
            nFft = obj.Nfft;
            
            % Standard window and its time derivative dw/dt
            w = chirp.dsp.WindowFunction.create('hann', winLen);
            tWin = (-(winLen - 1)/2 : (winLen - 1)/2).' / fs;
            % Analytical derivative of Hann window
            dw = -pi * sin(2.0 * pi * (0:winLen-1).' / (winLen - 1)) / ((winLen - 1) / fs);
            
            numFrames = floor((N - winLen) / hop) + 1;
            sstMatrix = zeros(nFft, numFrames);
            
            for k = 1:numFrames
                startIdx = (k - 1) * hop + 1;
                seg = x(startIdx:(startIdx + winLen - 1)).';
                
                % STFT with window w and derivative window dw
                Vg = fft(seg .* w, nFft);
                Vdg = fft(seg .* dw, nFft);
                
                % Instantaneous frequency estimate: omega_inst = real( Vdg / (2*pi*j * Vg) )
                omegaInst = imag(Vdg ./ max(abs(Vg), 1e-12));
                
                % Reallocate energy into sstMatrix bins
                for fBin = 1:nFft
                    val = abs(Vg(fBin));
                    if val > 1e-6
                        targetBin = mod(round(fBin + omegaInst(fBin) / (fs / nFft)) - 1, nFft) + 1;
                        sstMatrix(targetBin, k) = sstMatrix(targetBin, k) + val;
                    end
                end
            end
            
            sstMatrix = fftshift(sstMatrix, 1);
            timeAxis = ((0:numFrames-1) * hop + winLen / 2.0) / fs;
            
            if mod(nFft, 2) == 0
                kIndices = (-nFft/2):(nFft/2 - 1);
            else
                kIndices = (-(nFft - 1)/2):((nFft - 1)/2);
            end
            freqAxis = kIndices * (fs / nFft);
        end
    end
end
