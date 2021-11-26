classdef PhaseCodeModulator
    % PHASECODEMODULATOR Phase-coded waveforms for MIMO radar and spread-spectrum.
    %
    % Generates Barker sequences, Frank polyphase codes, and Zadoff-Chu CAZAC
    % sequences for orthogonal waveform synthesis and pulse compression.
    
    methods (Static)
        function code = barker(lengthCode)
            % Generate binary Barker sequence of specified length: 2, 3, 4, 5, 7, 11, or 13.
            % Returns a bipolar sequence with values +1 and -1.
            
            switch lengthCode
                case 2
                    code = [1, -1];
                case 3
                    code = [1, 1, -1];
                case 4
                    code = [1, 1, -1, 1];
                case 5
                    code = [1, 1, 1, -1, 1];
                case 7
                    code = [1, 1, 1, -1, -1, 1, -1];
                case 11
                    code = [1, 1, 1, -1, -1, -1, 1, -1, -1, 1, -1];
                case 13
                    code = [1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1];
                otherwise
                    error('PhaseCodeModulator:InvalidBarkerLength', ...
                        'Supported Barker code lengths are 2, 3, 4, 5, 7, 11, and 13.');
            end
        end
        
        function code = frank(M)
            % Generate Frank polyphase sequence of length M^2.
            %
            % Phase values: phi_{p, q} = (2*pi / M) * (p - 1) * (q - 1)
            % for p, q in 1:M.
            
            if M <= 1 || round(M) ~= M
                error('PhaseCodeModulator:InvalidFrankOrder', 'Frank order M must be an integer >= 2.');
            end
            
            code = zeros(1, M * M);
            idx = 1;
            for p = 1:M
                for q = 1:M
                    phi = (2.0 * pi / M) * (p - 1) * (q - 1);
                    code(idx) = exp(1j * phi);
                    idx = idx + 1;
                end
            end
        end
        
        function code = zadoff_chu(u, N_zc)
            % Generate Zadoff-Chu CAZAC complex sequence.
            %
            % Inputs:
            %   u    - Root index (must be coprime to N_zc)
            %   N_zc - Sequence length (prime length recommended)
            %
            % Sequence formula:
            %   x[n] = exp(-1j * pi * u * n * (n + mod(N_zc, 2)) / N_zc)
            
            if gcd(u, N_zc) ~= 1
                warning('PhaseCodeModulator:NotCoprime', 'Root u and length N_zc should be coprime for ideal CAZAC properties.');
            end
            
            n = 0:(N_zc - 1);
            q = mod(N_zc, 2);
            phase = -pi * u * n .* (n + q) / N_zc;
            code = exp(1j * phase);
        end
        
        function [r, lags] = aperiodic_autocorr(code)
            % Compute aperiodic discrete autocorrelation of code.
            code = code(:).';
            N = length(code);
            r = conv(code, fliplr(conj(code)));
            lags = -(N - 1):(N - 1);
        end
        
        function psl = peak_sidelobe_ratio_db(code)
            % Compute Peak Sidelobe Level (PSL) in dB.
            r = chirp.waveform.PhaseCodeModulator.aperiodic_autocorr(code);
            centerIdx = ceil(length(r) / 2);
            peakVal = abs(r(centerIdx));
            
            % Zero out main peak to isolate sidelobes
            r(centerIdx) = 0.0;
            maxSidelobe = max(abs(r));
            
            if maxSidelobe <= 0
                psl = -Inf;
            else
                psl = 20.0 * log10(maxSidelobe / peakVal);
            end
        end
    end
end
