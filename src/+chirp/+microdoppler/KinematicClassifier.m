classdef KinematicClassifier
    % KINEMATICCLASSIFIER Radar target classification based on micro-Doppler signatures.
    %
    % Classifies physical moving obstacles into target categories:
    %   'pedestrian', 'cyclist', 'vehicle', 'drone', or 'stationary_clutter'
    % using extracted bandwidth, cadence, entropy, and translational velocity features.
    
    methods (Static)
        function [className, confidence, classScores] = classify(feats, wavelengthMeters)
            % Classify target from extracted MicroDopplerFeatures struct.
            %
            % Inputs:
            %   feats            - Struct returned by chirp.microdoppler.MicroDopplerFeatures
            %   wavelengthMeters - Radar carrier wavelength (default: 3.89 mm for 77 GHz)
            %
            % Returns:
            %   className        - String category: 'pedestrian', 'cyclist', 'vehicle', 'drone'
            %   confidence       - Classification confidence score in [0, 1]
            %   classScores      - Struct with likelihood scores for all classes
            
            if nargin < 2 || isempty(wavelengthMeters)
                wavelengthMeters = chirp.utils.PhysicalConstants.freq2wavelength(77e9);
            end
            
            % Convert Doppler frequencies to radial velocities: v = f * lambda / 2
            vBase = abs(feats.torsoVelocityHz) * wavelengthMeters / 2.0;
            vBw = feats.meanBandwidthHz * wavelengthMeters / 2.0;
            cadence = feats.cadenceHz;
            entropy = feats.spectralEntropy;
            
            % 1. Vehicle: High torso velocity (> 6 m/s), low micro-Doppler bandwidth (< 1.0 m/s, rigid body)
            scoreVehicle = -((vBase - 15.0)^2) / (2 * 8.0^2) - ((vBw - 0.4)^2) / (2 * 0.4^2);
            
            % 2. Pedestrian: Slow translational velocity (0.5..2.5 m/s), limb bandwidth (0.8..2.5 m/s)
            scorePedestrian = -((vBase - 1.4)^2) / (2 * 0.8^2) - ((vBw - 1.5)^2) / (2 * 0.8^2);
            if cadence > 0.5 && cadence < 4.0
                scorePedestrian = scorePedestrian - ((cadence - 1.8)^2) / (2 * 0.8^2);
            end
            
            % 3. Cyclist: Medium velocity (3..8 m/s), pedaling cadence, wheel micro-Doppler (2..4 m/s)
            scoreCyclist = -((vBase - 5.5)^2) / (2 * 1.8^2) - ((vBw - 3.2)^2) / (2 * 1.2^2);
            if cadence > 0.5 && cadence < 4.0
                scoreCyclist = scoreCyclist - ((cadence - 2.2)^2) / (2 * 0.8^2);
            end
            
            % 4. Drone / UAV: High rotor blade tip speed (vBw > 8 m/s), high spectral entropy
            scoreDrone = -((vBw - 12.0)^2) / (2 * 4.0^2) - ((entropy - 0.90)^2) / (2 * 0.1^2);
            
            rawScores = [scoreVehicle, scorePedestrian, scoreCyclist, scoreDrone];
            maxScore = max(rawScores);
            expScores = exp(rawScores - maxScore);
            probs = expScores / sum(expScores);
            
            classes = {'vehicle', 'pedestrian', 'cyclist', 'drone'};
            [conf, bestIdx] = max(probs);
            className = classes{bestIdx};
            confidence = conf;
            
            classScores = struct(...
                'vehicle', probs(1), ...
                'pedestrian', probs(2), ...
                'cyclist', probs(3), ...
                'drone', probs(4) ...
            );
        end
    end
end
