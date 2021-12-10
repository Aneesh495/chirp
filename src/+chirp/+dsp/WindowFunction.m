classdef WindowFunction
    % WINDOWFUNCTION Spectral window generation and noise equivalent bandwidth calculator.
    %
    % Implements Hann, Hamming, Blackman, Blackman-Harris, Flat-Top, Kaiser,
    % and Chebyshev windows from first principles without external toolbox dependencies.
    
    methods (Static)
        function w = create(windowType, N, param)
            % Generate an N-point window vector.
            %
            % Inputs:
            %   windowType - String: 'rect', 'hann', 'hamming', 'blackman',
            %                'blackman-harris', 'flattop', 'kaiser', 'chebyshev'
            %   N          - Integer window length
            %   param      - Optional tuning parameter (e.g. beta for Kaiser, r for Chebyshev)
            %
            % Output:
            %   w          - N x 1 column vector of window weights
            
            if N <= 0
                error('WindowFunction:InvalidLength', 'Window length must be strictly positive.');
            end
            
            if N == 1
                w = 1.0;
                return;
            end
            
            n = (0:(N - 1)).';
            
            switch lower(windowType)
                case 'rect'
                    w = ones(N, 1);
                    
                case 'hann'
                    w = 0.5 - 0.5 * cos(2.0 * pi * n / (N - 1));
                    
                case 'hamming'
                    w = 0.54 - 0.46 * cos(2.0 * pi * n / (N - 1));
                    
                case 'blackman'
                    w = 0.42 - 0.5 * cos(2.0 * pi * n / (N - 1)) + 0.08 * cos(4.0 * pi * n / (N - 1));
                    
                case 'blackman-harris'
                    % 4-term Blackman-Harris window (-92 dB sidelobe suppression)
                    a0 = 0.35875;
                    a1 = 0.48829;
                    a2 = 0.14128;
                    a3 = 0.01168;
                    w = a0 - a1 * cos(2.0 * pi * n / (N - 1)) + ...
                             a2 * cos(4.0 * pi * n / (N - 1)) - ...
                             a3 * cos(6.0 * pi * n / (N - 1));
                             
                case 'flattop'
                    % Flat-top window for high-precision amplitude calibration
                    a0 = 0.21557895;
                    a1 = 0.41663158;
                    a2 = 0.277263158;
                    a3 = 0.083578947;
                    a4 = 0.006947368;
                    w = a0 - a1 * cos(2.0 * pi * n / (N - 1)) + ...
                             a2 * cos(4.0 * pi * n / (N - 1)) - ...
                             a3 * cos(6.0 * pi * n / (N - 1)) + ...
                             a4 * cos(8.0 * pi * n / (N - 1));
                             
                case 'kaiser'
                    if nargin < 3 || isempty(param), beta = 7.865; else, beta = param; end
                    alpha = (N - 1) / 2.0;
                    r = (n - alpha) / alpha;
                    arg = beta * sqrt(max(1.0 - r.^2, 0.0));
                    w = besseli(0, arg) ./ besseli(0, beta);
                    
                case 'chebyshev'
                    % Dolph-Chebyshev window with specified sidelobe attenuation in dB
                    if nargin < 3 || isempty(param), rDb = 80.0; else, rDb = param; end
                    w = chirp.dsp.WindowFunction.dolphChebyshev(N, rDb);
                    
                otherwise
                    error('WindowFunction:UnknownType', 'Window type %s not supported.', windowType);
            end
        end
        
        function metrics = evaluateMetrics(w)
            % Compute spectral figures of merit for window vector w.
            %
            % Returns:
            %   metrics.coherentGain - Sum(w) / N
            %   metrics.enbw         - N * sum(w.^2) / (sum(w)^2) (Equivalent Noise Bandwidth in bins)
            %   metrics.scallopingDb - Scalloping loss at bin boundary (dB)
            
            w = w(:);
            N = length(w);
            sumW = sum(w);
            sumW2 = sum(w.^2);
            
            coherentGain = sumW / N;
            enbw = N * sumW2 / (sumW^2);
            
            % Half-bin frequency response for scalloping loss
            halfBinPhase = exp(-1j * pi * (0:(N - 1)).' / N);
            valHalfBin = abs(sum(w .* halfBinPhase)) / sumW;
            scallopingDb = -20.0 * log10(max(valHalfBin, 1e-12));
            
            metrics = struct(...
                'coherentGain', coherentGain, ...
                'enbw', enbw, ...
                'scallopingDb', scallopingDb ...
            );
        end
    end
    
    methods (Static, Access = private)
        function w = dolphChebyshev(N, rDb)
            % Exact Dolph-Chebyshev window synthesis via Chebyshev polynomials.
            M = N - 1;
            r = 10^(rDb / 20.0);
            x0 = cosh(acosh(r) / M);
            
            k = 0:M;
            theta = pi * k / N;
            x = x0 * cos(theta);
            
            % Evaluate Chebyshev polynomial T_M(x)
            W = zeros(1, N);
            for iK = 1:N
                xi = x(iK);
                if abs(xi) <= 1.0
                    W(iK) = cos(M * acos(xi));
                else
                    W(iK) = cosh(M * acosh(xi));
                end
            end
            
            % Inverse DFT to obtain time-domain window weights
            w = real(ifft(W));
            w = fftshift(w);
            w = w(:) / max(w);
        end
    end
end
