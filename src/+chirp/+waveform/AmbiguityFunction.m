classdef AmbiguityFunction
    % AMBIGUITYFUNCTION Woodward ambiguity function evaluator for radar waveforms.
    %
    % Computes the 2D Woodward ambiguity surface |chi(tau, nu)|^2 and 1D cuts
    % to analyze range resolution, Doppler tolerance, and coupling sidelobes.
    
    methods (Static)
        function [af, tauGrid, nuGrid] = compute(signal, fs, maxDelay, maxDoppler, numTau, numNu)
            % Compute discrete 2D Woodward ambiguity function for a finite complex signal.
            %
            % Inputs:
            %   signal     - 1D complex discrete signal vector
            %   fs         - Sampling rate (Hz)
            %   maxDelay   - Maximum delay tau to evaluate (s)
            %   maxDoppler - Maximum Doppler frequency nu to evaluate (Hz)
            %   numTau     - Number of delay evaluation points
            %   numNu      - Number of Doppler evaluation points
            %
            % Returns:
            %   af         - numNu x numTau normalized ambiguity magnitude matrix
            %   tauGrid    - 1 x numTau delay evaluation vector (s)
            %   nuGrid     - 1 x numNu Doppler evaluation vector (Hz)
            
            signal = signal(:).';
            N = length(signal);
            
            if nargin < 3 || isempty(maxDelay), maxDelay = (N / fs) * 0.5; end
            if nargin < 4 || isempty(maxDoppler), maxDoppler = fs * 0.25; end
            if nargin < 5 || isempty(numTau), numTau = 101; end
            if nargin < 6 || isempty(numNu), numNu = 101; end
            
            tauGrid = linspace(-maxDelay, maxDelay, numTau);
            nuGrid = linspace(-maxDoppler, maxDoppler, numNu);
            
            af = zeros(numNu, numTau);
            t = (0:(N - 1)) ./ fs;
            
            % Compute energy of the signal for normalization
            eSig = sum(abs(signal).^2);
            if eSig <= 0, eSig = 1.0; end
            
            for iTau = 1:numTau
                tau = tauGrid(iTau);
                delaySamples = round(tau * fs);
                
                % Form shifted signal s(t - tau)
                if delaySamples >= 0
                    sigShifted = [zeros(1, delaySamples), signal(1:(N - delaySamples))];
                else
                    sigShifted = [signal((abs(delaySamples) + 1):N), zeros(1, abs(delaySamples))];
                end
                
                % Cross-product term s(t) * conj(s(t - tau))
                prodSig = signal .* conj(sigShifted);
                
                for iNu = 1:numNu
                    nu = nuGrid(iNu);
                    carrier = exp(1j * 2.0 * pi * nu * t);
                    val = sum(prodSig .* carrier);
                    af(iNu, iTau) = abs(val);
                end
            end
            
            % Normalize peak to 1.0 (0 dB)
            af = af ./ eSig;
        end
        
        function [delayCut, tauGrid] = zeroDopplerCut(signal, fs, maxDelay, numTau)
            % Compute 1D zero-Doppler cut (autocorrelation function).
            if nargin < 3, maxDelay = []; end
            if nargin < 4, numTau = 201; end
            [af, tauGrid, ~] = chirp.waveform.AmbiguityFunction.compute(signal, fs, maxDelay, 0, numTau, 1);
            delayCut = af(1, :);
        end
        
        function [dopplerCut, nuGrid] = zeroDelayCut(signal, fs, maxDoppler, numNu)
            % Compute 1D zero-delay cut (spectral energy magnitude).
            if nargin < 3, maxDoppler = []; end
            if nargin < 4, numNu = 201; end
            [af, ~, nuGrid] = chirp.waveform.AmbiguityFunction.compute(signal, fs, 0, maxDoppler, 1, numNu);
            dopplerCut = af(:, 1).';
        end
    end
end
