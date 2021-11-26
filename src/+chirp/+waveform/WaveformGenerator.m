classdef WaveformGenerator
    % WAVEFORMGENERATOR Factory for standardized FMCW and NLFM radar waveforms.
    %
    % Creates standard automotive 77 GHz LRR, MRR, USRR profiles, airborne
    % SAR chirp configurations, and tangent-based non-linear chirp waveforms.
    
    methods (Static)
        function chirpObj = automotiveLRR()
            % Automotive Long-Range Radar (LRR) chirp profile:
            %   fc = 77.0 GHz
            %   B  = 500 MHz (deltaR = 0.30 m)
            %   Tc = 50.0 us
            %   Fs = 20.0 MHz
            %   Rmax = 150 m
            chirpObj = chirp.waveform.FmcwChirp(...
                'StartFrequency', 77.0e9, ...
                'SweepBandwidth', 500.0e6, ...
                'ChirpDuration', 50.0e-6, ...
                'SampleRate', 20.0e6, ...
                'IdleTime', 7.0e-6, ...
                'AdcStartTime', 3.0e-6, ...
                'NumAdcSamples', 512 ...
            );
        end
        
        function chirpObj = automotiveSRR()
            % Automotive Short-Range Radar (SRR) / High-Resolution Parking chirp profile:
            %   fc = 77.0 GHz
            %   B  = 2.0 GHz (deltaR = 0.075 m)
            %   Tc = 30.0 us
            %   Fs = 30.0 MHz
            chirpObj = chirp.waveform.FmcwChirp(...
                'StartFrequency', 77.0e9, ...
                'SweepBandwidth', 2.0e9, ...
                'ChirpDuration', 30.0e-6, ...
                'SampleRate', 30.0e6, ...
                'IdleTime', 5.0e-6, ...
                'AdcStartTime', 2.0e-6, ...
                'NumAdcSamples', 512 ...
            );
        end
        
        function chirpObj = airborneSAR()
            % Airborne X-Band Synthetic Aperture Radar (SAR) chirp profile:
            %   fc = 9.6 GHz
            %   B  = 600 MHz
            %   Tc = 20.0 us
            %   Fs = 50.0 MHz
            chirpObj = chirp.waveform.FmcwChirp(...
                'StartFrequency', 9.6e9, ...
                'SweepBandwidth', 600.0e6, ...
                'ChirpDuration', 20.0e-6, ...
                'SampleRate', 50.0e6, ...
                'IdleTime', 10.0e-6, ...
                'AdcStartTime', 1.0e-6, ...
                'NumAdcSamples', 1024 ...
            );
        end
        
        function [sig, t] = nlfmTangent(B, Tc, Fs, alphaParam)
            % Generate Non-Linear Frequency Modulated (NLFM) tangent chirp waveform.
            %
            % Modulates instantaneous frequency using a tangent trajectory:
            %   f(t) = (B / (2 * tan(alpha))) * tan(2 * alpha * t / Tc)
            % for t in [-Tc/2, Tc/2].
            %
            % Sinks more energy near edges of band to naturally suppress range sidelobes
            % without requiring window weighting that lowers SNR.
            
            if nargin < 4 || isempty(alphaParam), alphaParam = 1.35; end
            
            N = round(Tc * Fs);
            t = linspace(-Tc / 2.0, Tc / 2.0, N);
            
            % Frequency trajectory
            normT = (2.0 * alphaParam / Tc) * t;
            freq = (B / (2.0 * tan(alphaParam))) * tan(normT);
            
            % Integrate frequency numerically to obtain instantaneous phase
            dt = t(2) - t(1);
            phase = 2.0 * pi * cumsum(freq) * dt;
            
            sig = exp(1j * phase);
        end
    end
end
