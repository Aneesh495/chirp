classdef WignerVilleDistribution
    % WIGNERVILLEDISTRIBUTION High-resolution quadratic time-frequency distribution.
    %
    % Computes the Smoothed Pseudo-Wigner-Ville Distribution (SPWVD) with separable
    % time and frequency smoothing windows to eliminate non-linear cross-terms.
    
    properties
        TimeWindowLength = 31;        % Length of time smoothing window g(t)
        FreqWindowLength = 63;        % Length of frequency smoothing window h(tau)
        Nfft = 128;                   % FFT evaluation length
    end
    
    methods
        function obj = WignerVilleDistribution(varargin)
            % Construct a WignerVilleDistribution instance with optional parameters.
            if nargin > 0
                for k = 1:2:nargin
                    propName = varargin{k};
                    propVal = varargin{k+1};
                    if isprop(obj, propName)
                        obj.(propName) = propVal;
                    else
                        error('WignerVilleDistribution:UnknownProperty', 'Property %s not found.', propName);
                    end
                end
            end
        end
        
        function [spwvdDb, timeAxis, freqAxis] = compute(obj, signal, fs)
            % Compute Smoothed Pseudo-Wigner-Ville Distribution.
            %
            % Formula:
            %   W(t, f) = sum_tau [ h(tau) * sum_u (g(u) * x(t+u+tau) * conj(x(t+u-tau))) ] * e^(-j*4*pi*f*tau)
            
            x = signal(:).';
            N = length(x);
            nFft = obj.Nfft;
            
            % Smoothing windows
            g = chirp.dsp.WindowFunction.create('hamming', obj.TimeWindowLength);
            g = (g / sum(g)).';
            h = chirp.dsp.WindowFunction.create('hamming', obj.FreqWindowLength);
            h = (h / sum(h)).';
            
            halfG = floor(obj.TimeWindowLength / 2);
            halfH = floor(obj.FreqWindowLength / 2);
            
            W = zeros(nFft, N);
            
            for t = (halfG + halfH + 1):(N - halfG - halfH)
                tauVec = zeros(1, nFft);
                
                for tau = -halfH:halfH
                    % Time smoothed instantaneous autocorrelation R_t(tau)
                    r_val = 0.0;
                    for u = -halfG:halfG
                        idx1 = t + u + tau;
                        idx2 = t + u - tau;
                        r_val = r_val + g(u + halfG + 1) * x(idx1) * conj(x(idx2));
                    end
                    
                    tauIdx = mod(tau, nFft) + 1;
                    tauVec(tauIdx) = h(tau + halfH + 1) * r_val;
                end
                
                col = real(fftshift(fft(tauVec, nFft)));
                W(:, t) = col(:);
            end
            
            spwvdDb = 10.0 * log10(max(W, eps));
            timeAxis = (0:(N - 1)) / fs;
            
            if mod(nFft, 2) == 0
                kIndices = (-nFft/2):(nFft/2 - 1);
            else
                kIndices = (-(nFft - 1)/2):((nFft - 1)/2);
            end
            freqAxis = kIndices * (fs / nFft);
        end
    end
end
