classdef FmcwChirp
    % FMCWCHIRP Linear frequency modulated continuous wave (FMCW) chirp model.
    %
    % Models the physical transmit chirp parameters, fast-time sampling grid,
    % analytical baseband signal, and de-chirped beat signal formation.
    
    properties
        StartFrequency = 77.0e9;       % Carrier frequency (Hz)
        SweepBandwidth = 1.0e9;        % Sweep bandwidth B (Hz)
        ChirpDuration = 40.0e-6;       % Chirp active ramp duration Tc (s)
        SampleRate = 25.0e6;           % ADC sampling rate Fs (Hz)
        IdleTime = 5.0e-6;             % Inter-chirp idle time (s)
        AdcStartTime = 2.0e-6;         % ADC sampling delay after ramp start (s)
        NumAdcSamples = 512;           % Fast-time sample count
    end
    
    properties (Dependent)
        Slope                 % Chirp frequency slope S = B / Tc (Hz/s)
        Wavelength            % Carrier wavelength lambda = c / fc (m)
        RangeResolution       % Range resolution deltaR = c / (2 * B) (m)
        MaxUnambiguousRange   % Max beat range Rmax = (c * Fs) / (4 * S) (m)
        TotalChirpPeriod      % Tchirp = ChirpDuration + IdleTime (s)
    end
    
    methods
        function obj = FmcwChirp(varargin)
            % Construct an FMCW chirp configuration with optional name-value pairs.
            if nargin > 0
                for k = 1:2:nargin
                    propName = varargin{k};
                    propVal = varargin{k+1};
                    if isprop(obj, propName)
                        obj.(propName) = propVal;
                    else
                        error('FmcwChirp:UnknownProperty', 'Property %s not found in FmcwChirp.', propName);
                    end
                end
            end
        end
        
        function s = get.Slope(obj)
            s = obj.SweepBandwidth / obj.ChirpDuration;
        end
        
        function w = get.Wavelength(obj)
            w = chirp.utils.PhysicalConstants.C / obj.StartFrequency;
        end
        
        function dr = get.RangeResolution(obj)
            dr = chirp.utils.PhysicalConstants.C / (2.0 * obj.SweepBandwidth);
        end
        
        function rmax = get.MaxUnambiguousRange(obj)
            % Beat frequency fb = (2 * S * R) / c <= Fs / 2
            % Therefore Rmax = (c * Fs) / (4 * S)
            rmax = (chirp.utils.PhysicalConstants.C * obj.SampleRate) / (4.0 * obj.Slope);
        end
        
        function tp = get.TotalChirpPeriod(obj)
            tp = obj.ChirpDuration + obj.IdleTime;
        end
        
        function t = fastTimeVector(obj)
            % Returns fast-time sampling vector within ADC acquisition window.
            t = obj.AdcStartTime + (0:(obj.NumAdcSamples - 1)) ./ obj.SampleRate;
        end
        
        function [beatSignal, fb] = dechirpTarget(obj, rangeMeters, velocityMps, amplitude)
            % Compute de-chirped IF beat signal for a point target at range R and radial velocity v.
            %
            % Inputs:
            %   rangeMeters - Radial range to target (m)
            %   velocityMps - Radial target velocity along line-of-sight (m/s, positive = receding)
            %   amplitude   - Complex received amplitude scalar
            %
            % Returns:
            %   beatSignal  - 1 x NumAdcSamples complex double vector
            %   fb          - Analytical beat frequency (Hz)
            
            if nargin < 4 || isempty(amplitude)
                amplitude = 1.0;
            end
            
            c = chirp.utils.PhysicalConstants.C;
            tau = 2.0 * rangeMeters / c;  % Round-trip propagation delay
            
            % Beat frequency fb = (2 * S * R)/c + (2 * fc * v)/c
            fb_range = obj.Slope * tau;
            fb_doppler = (2.0 * obj.StartFrequency * velocityMps) / c;
            fb = fb_range + fb_doppler;
            
            t = obj.fastTimeVector();
            
            % Exact de-chirped phase:
            % phi_b(t) = 2*pi * (fb * t + fc * tau - 0.5 * S * tau^2)
            phi_b = 2.0 * pi * (fb * t + obj.StartFrequency * tau - 0.5 * obj.Slope * (tau^2));
            beatSignal = amplitude * exp(1j * phi_b);
        end
    end
end
