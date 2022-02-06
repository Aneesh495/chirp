classdef CadenceExtractor
    % CADENCEEXTRACTOR Periodicity and gait cadence estimator for micro-Doppler signals.
    %
    % Uses cepstral peak detection and autocorrelation methods to measure stride
    % frequency in pedestrians, pedal cadence in cyclists, and rotor RPM in drones.
    
    methods (Static)
        function [cadenceHz, confidence] = extractFromEnvelope(envelopeSignal, fs)
            % Estimate fundamental periodicity from envelope time-series.
            %
            % Inputs:
            %   envelopeSignal - 1D real vector of time-varying Doppler extent
            %   fs             - Sampling rate of the envelope (Hz)
            %
            % Returns:
            %   cadenceHz      - Estimated dominant frequency in Hz
            %   confidence     - Periodic peak prominence metric in [0, 1]
            
            x = envelopeSignal(:).';
            N = length(x);
            if N < 16
                cadenceHz = 0.0;
                confidence = 0.0;
                return;
            end
            
            % Remove DC trend
            x = x - mean(x);
            
            % Autocorrelation
            r = conv(x, fliplr(x));
            mid = ceil(length(r) / 2);
            rHalf = r(mid:end);
            rNorm = rHalf / max(rHalf(1), eps);
            
            % Search within plausible biological/mechanical cadence window [0.5 Hz, 8.0 Hz]
            minLag = max(2, round(fs / 8.0));
            maxLag = min(length(rNorm) - 1, round(fs / 0.5));
            
            if minLag >= maxLag
                cadenceHz = 0.0;
                confidence = 0.0;
                return;
            end
            
            searchRegion = rNorm(minLag:maxLag);
            [pkVal, relIdx] = max(searchRegion);
            bestLag = minLag + relIdx - 1;
            
            cadenceHz = fs / bestLag;
            confidence = max(0.0, min(1.0, pkVal));
        end
    end
end
