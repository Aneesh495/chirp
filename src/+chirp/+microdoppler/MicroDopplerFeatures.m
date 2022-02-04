classdef MicroDopplerFeatures
    % MICRODOPPLERFEATURES Biomechanical and kinematic feature extraction from spectrograms.
    %
    % Extracts torso velocity envelope, micro-Doppler bandwidth, cadence periodicity,
    % spectral entropy, and blade flash signatures from time-frequency energy distributions.
    
    methods (Static)
        function feats = extract(spectrogramPower, timeAxis, freqAxis)
            % Extract scalar and profile kinematic features from a 2D spectrogram.
            %
            % Inputs:
            %   spectrogramPower - [N_freq x N_time] linear power matrix
            %   timeAxis         - 1 x N_time vector (s)
            %   freqAxis         - 1 x N_freq vector (Hz)
            %
            % Returns:
            %   feats struct containing:
            %     meanBandwidthHz  - Average micro-Doppler bandwidth extent
            %     maxDopplerHz     - Peak Doppler excursion
            %     cadenceHz        - Dominant gait/pedaling/rotation cadence
            %     torsoVelocityHz  - Central energy envelope
            %     spectralEntropy  - Normalized Shannon entropy across Doppler bands
            
            [nFreq, nTime] = size(spectrogramPower);
            
            % Compute upper and lower frequency envelopes at threshold (-15 dB from column peak)
            upperEnv = zeros(1, nTime);
            lowerEnv = zeros(1, nTime);
            centerEnv = zeros(1, nTime);
            
            for t = 1:nTime
                colPwr = spectrogramPower(:, t);
                pkVal = max(colPwr);
                thresh = pkVal * 0.05; % -13 dB cutoff
                
                sigBins = find(colPwr >= thresh);
                if isempty(sigBins)
                    upperEnv(t) = 0.0;
                    lowerEnv(t) = 0.0;
                    centerEnv(t) = 0.0;
                else
                    upperEnv(t) = max(freqAxis(sigBins));
                    lowerEnv(t) = min(freqAxis(sigBins));
                    % Energy centroid (spectral mean)
                    centerEnv(t) = sum(freqAxis(sigBins).' .* colPwr(sigBins)) / max(sum(colPwr(sigBins)), eps);
                end
            end
            
            bandwidth = upperEnv - lowerEnv;
            meanBw = mean(bandwidth);
            
            [pkRow, ~] = find(spectrogramPower == max(spectrogramPower(:)), 1);
            maxDop = abs(freqAxis(pkRow));
            
            % Cadence extraction via autocorrelation of bandwidth envelope
            bwZeroMean = bandwidth - mean(bandwidth);
            if std(bwZeroMean) > 1e-6
                rBw = conv(bwZeroMean, fliplr(bwZeroMean));
                centerIdx = ceil(length(rBw) / 2);
                halfR = rBw(centerIdx:end);
                
                % Find first dominant peak away from zero-lag
                diffSigns = diff(sign(diff(halfR)));
                peakLags = find(diffSigns < 0);
                if ~isempty(peakLags)
                    bestLag = peakLags(1) + 1;
                    dt = timeAxis(2) - timeAxis(1);
                    periodSec = bestLag * dt;
                    cadence = 1.0 / max(periodSec, 1e-3);
                else
                    cadence = 0.0;
                end
            else
                cadence = 0.0;
            end
            
            % Spectral entropy: measure of energy dispersion across Doppler bins
            pNorm = spectrogramPower / max(sum(spectrogramPower(:)), eps);
            pNorm = max(pNorm, eps);
            entropyVal = -sum(pNorm(:) .* log(pNorm(:))) / log(nFreq * nTime);
            
            feats = struct(...
                'meanBandwidthHz', meanBw, ...
                'maxDopplerHz', maxDop, ...
                'cadenceHz', cadence, ...
                'torsoVelocityHz', mean(centerEnv), ...
                'spectralEntropy', entropyVal, ...
                'upperEnvelope', upperEnv, ...
                'lowerEnvelope', lowerEnv ...
            );
        end
    end
end
