classdef TmCfar
    % TMCFAR Trimmed-Mean Constant False Alarm Rate detector.
    %
    % Sorts reference training cells and discards T1 smallest and T2 largest values,
    % averaging the central subset. Balances low CFAR SNR loss with target masking immunity.
    
    properties
        Pfa = 1e-4;                   % False alarm probability
        NumTrain = 20;                % Total training cells (sum of both sides)
        NumGuard = 4;                 % Total guard cells
        TrimLow = 2;                  % Number of smallest values to discard
        TrimHigh = 2;                 % Number of largest values to discard
    end
    
    methods
        function obj = TmCfar(varargin)
            % Construct a TmCfar detector instance.
            if nargin > 0
                for k = 1:2:nargin
                    propName = varargin{k};
                    propVal = varargin{k+1};
                    if isprop(obj, propName)
                        obj.(propName) = propVal;
                    else
                        error('TmCfar:UnknownProperty', 'Property %s not recognized.', propName);
                    end
                end
            end
        end
        
        function [detectionMask, thresholdVec] = detect1D(obj, signalPower)
            % 1D Trimmed-Mean CFAR detection.
            x = signalPower(:).';
            N = length(x);
            nT = obj.NumTrain;
            nG = obj.NumGuard;
            
            nRemaining = nT - (obj.TrimLow + obj.TrimHigh);
            if nRemaining <= 0
                error('TmCfar:InvalidTrimming', 'Trim parameters leave zero reference cells.');
            end
            
            % Nominal scaling factor derived from remaining training cells
            alpha = chirp.cfar.ThresholdCalculator.caCfarAlpha(obj.Pfa, nRemaining);
            
            thresholdVec = zeros(1, N);
            halfWin = round((nT + nG) / 2);
            halfG = round(nG / 2);
            
            for i = 1:N
                leftIdx = max(1, i - halfWin) : max(1, i - halfG - 1);
                rightIdx = min(N, i + halfG + 1) : min(N, i + halfWin);
                
                combined = [x(leftIdx), x(rightIdx)];
                if length(combined) > (obj.TrimLow + obj.TrimHigh)
                    sortedVals = sort(combined);
                    validSubset = sortedVals((obj.TrimLow + 1):(end - obj.TrimHigh));
                    nf = mean(validSubset);
                else
                    nf = mean(combined);
                end
                
                thresholdVec(i) = alpha * max(nf, 1e-12);
            end
            
            detectionMask = (x > thresholdVec);
        end
    end
end
