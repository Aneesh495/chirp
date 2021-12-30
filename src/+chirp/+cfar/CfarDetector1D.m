classdef CfarDetector1D
    % CFARDETECTOR1D 1D adaptive Constant False Alarm Rate sliding window detector.
    %
    % Performs adaptive noise floor estimation and target thresholding along 1D
    % range profiles or Doppler slices. Supports CA, OS, GO, and SO modes.
    
    properties
        Mode = 'CA';                  % Detection algorithm: 'CA', 'OS', 'GO', 'SO'
        Pfa = 1e-4;                   % Desired probability of false alarm
        NumTrain = 16;                % Number of one-sided training cells
        NumGuard = 2;                 % Number of one-sided guard cells
        RankK = [];                   % OS-CFAR rank statistic index
    end
    
    methods
        function obj = CfarDetector1D(varargin)
            % Construct CfarDetector1D instance with optional configuration.
            if nargin > 0
                for k = 1:2:nargin
                    propName = varargin{k};
                    propVal = varargin{k+1};
                    if isprop(obj, propName)
                        obj.(propName) = propVal;
                    else
                        error('CfarDetector1D:UnknownProperty', 'Property %s not found.', propName);
                    end
                end
            end
        end
        
        function [detectionMask, thresholdVector, noiseFloorVector] = detect(obj, signalPower)
            % Execute 1D CFAR detection across a linear power signal vector.
            %
            % Inputs:
            %   signalPower        - 1D vector of non-negative power values
            %
            % Returns:
            %   detectionMask      - Boolean vector (true where signalPower > threshold)
            %   thresholdVector    - Adaptive threshold values across all bins
            %   noiseFloorVector   - Estimated local clutter/noise floor power
            
            x = signalPower(:).';
            N = length(x);
            
            nT = obj.NumTrain;
            nG = obj.NumGuard;
            nTotalTrain = 2 * nT;
            
            % Compute scaling factor alpha
            switch upper(obj.Mode)
                case 'CA'
                    alpha = chirp.cfar.ThresholdCalculator.caCfarAlpha(obj.Pfa, nTotalTrain);
                case 'OS'
                    kRank = obj.RankK;
                    if isempty(kRank), kRank = round(0.75 * nTotalTrain); end
                    alpha = chirp.cfar.ThresholdCalculator.osCfarAlpha(obj.Pfa, nTotalTrain, kRank);
                case 'GO'
                    alpha = chirp.cfar.ThresholdCalculator.goCfarAlpha(obj.Pfa, nTotalTrain);
                case 'SO'
                    alpha = chirp.cfar.ThresholdCalculator.soCfarAlpha(obj.Pfa, nTotalTrain);
                otherwise
                    alpha = chirp.cfar.ThresholdCalculator.caCfarAlpha(obj.Pfa, nTotalTrain);
            end
            
            noiseFloorVector = zeros(1, N);
            thresholdVector = zeros(1, N);
            
            halfWin = nT + nG;
            
            for i = 1:N
                % Determine available window boundaries with boundary clamping
                leftTrainIdx = max(1, i - halfWin) : max(1, i - nG - 1);
                rightTrainIdx = min(N, i + nG + 1) : min(N, i + halfWin);
                
                % Filter out invalid / empty index ranges
                leftTrain = x(leftTrainIdx);
                rightTrain = x(rightTrainIdx);
                
                if isempty(leftTrain) && isempty(rightTrain)
                    noiseFloor = 1e-12;
                elseif isempty(leftTrain)
                    noiseFloor = mean(rightTrain);
                elseif isempty(rightTrain)
                    noiseFloor = mean(leftTrain);
                else
                    switch upper(obj.Mode)
                        case 'CA'
                            noiseFloor = (sum(leftTrain) + sum(rightTrain)) / (length(leftTrain) + length(rightTrain));
                        case 'OS'
                            combined = [leftTrain, rightTrain];
                            sortedVals = sort(combined);
                            kIdx = min(kRank, length(sortedVals));
                            noiseFloor = sortedVals(kIdx);
                        case 'GO'
                            noiseFloor = max(mean(leftTrain), mean(rightTrain));
                        case 'SO'
                            noiseFloor = min(mean(leftTrain), mean(rightTrain));
                        otherwise
                            noiseFloor = (sum(leftTrain) + sum(rightTrain)) / (length(leftTrain) + length(rightTrain));
                    end
                end
                
                noiseFloorVector(i) = noiseFloor;
                thresholdVector(i) = alpha * noiseFloor;
            end
            
            detectionMask = (x > thresholdVector);
        end
    end
end
