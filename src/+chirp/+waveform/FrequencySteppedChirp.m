classdef FrequencySteppedChirp
    % FREQUENCYSTEPPEDCHIRP Stepped-Frequency Continuous-Wave (SFCW) radar model.
    %
    % Emits a sequence of N pulses where carrier frequency steps by deltaF per pulse.
    % Reconstructs synthetic wideband high-resolution range profiles via IDFT.
    
    properties
        BaseFrequency = 76.0e9;        % Initial carrier frequency f0 (Hz)
        FrequencyStep = 10.0e6;        % Step size deltaF (Hz)
        NumSteps = 128;                % Number of steps N
        PulseDuration = 5.0e-6;        % Duration of each sub-pulse (s)
    end
    
    properties (Dependent)
        SyntheticBandwidth   % B_syn = N * deltaF (Hz)
        SyntheticRangeRes    % deltaR = c / (2 * B_syn) (m)
        MaxUnambiguousRange  % Rmax = c / (2 * deltaF) (m)
    end
    
    methods
        function obj = FrequencySteppedChirp(varargin)
            % Construct an SFCW configuration with optional name-value pairs.
            if nargin > 0
                for k = 1:2:nargin
                    propName = varargin{k};
                    propVal = varargin{k+1};
                    if isprop(obj, propName)
                        obj.(propName) = propVal;
                    else
                        error('FrequencySteppedChirp:UnknownProperty', 'Property %s not found.', propName);
                    end
                end
            end
        end
        
        function b = get.SyntheticBandwidth(obj)
            b = obj.NumSteps * obj.FrequencyStep;
        end
        
        function dr = get.SyntheticRangeRes(obj)
            dr = chirp.utils.PhysicalConstants.C / (2.0 * obj.SyntheticBandwidth);
        end
        
        function rmax = get.MaxUnambiguousRange(obj)
            rmax = chirp.utils.PhysicalConstants.C / (2.0 * obj.FrequencyStep);
        end
        
        function rxSamples = simulateReturns(obj, targetRanges, targetAmplitudes)
            % Simulate baseband I/Q samples collected across the N frequency steps.
            %
            % Inputs:
            %   targetRanges     - Vector of target distances (m)
            %   targetAmplitudes - Vector of complex target reflection amplitudes
            %
            % Returns:
            %   rxSamples        - 1 x NumSteps complex vector
            
            c = chirp.utils.PhysicalConstants.C;
            N = obj.NumSteps;
            rxSamples = zeros(1, N);
            
            stepIndices = 0:(N - 1);
            frequencies = obj.BaseFrequency + stepIndices * obj.FrequencyStep;
            
            for iTgt = 1:length(targetRanges)
                R = targetRanges(iTgt);
                amp = targetAmplitudes(iTgt);
                tau = 2.0 * R / c;
                % Round-trip phase shift: phi_n = 2*pi * f_n * tau
                phaseShift = exp(-1j * 2.0 * pi * frequencies * tau);
                rxSamples = rxSamples + amp * phaseShift;
            end
        end
        
        function [profile, rangeAxis] = reconstructProfile(obj, rxSamples, fftSize)
            % Reconstruct 1D range profile via Inverse Discrete Fourier Transform.
            %
            % Inputs:
            %   rxSamples - 1 x NumSteps received frequency-domain samples
            %   fftSize   - Zero-padded IDFT length (defaults to next power of 2)
            %
            % Returns:
            %   profile   - Complex range profile magnitude
            %   rangeAxis - Physical range coordinates (m)
            
            N = obj.NumSteps;
            if nargin < 3 || isempty(fftSize)
                fftSize = max(N, 2^chirp.utils.FastMath.nextpow2_fast(N * 4));
            end
            
            % Apply Hann window across frequency steps to reduce range sidelobes
            win = 0.5 - 0.5 * cos(2.0 * pi * (0:N-1) / (N - 1));
            windowedSamples = rxSamples(:).' .* win;
            
            profile = ifft(windowedSamples, fftSize);
            
            c = chirp.utils.PhysicalConstants.C;
            deltaR = c / (2.0 * obj.SyntheticBandwidth);
            rangeAxis = (0:(fftSize - 1)) * (obj.MaxUnambiguousRange / fftSize);
        end
    end
end
