classdef CfarDetector2D
    % CFARDETECTOR2D 2D adaptive Constant False Alarm Rate detector with integral image acceleration.
    %
    % Operates on 2D Range-Doppler matrices using rectangular training and guard windows.
    % Uses 2D cumulative sum integral tables for O(1) constant-time window evaluation.
    
    properties
        Mode = 'CA';                  % Detection algorithm: 'CA', 'OS', 'GO', 'SO'
        Pfa = 1e-4;                   % False alarm probability
        TrainRange = 8;               % One-sided training cells in range (rows)
        GuardRange = 2;               % One-sided guard cells in range (rows)
        TrainVelocity = 4;            % One-sided training cells in velocity (columns)
        GuardVelocity = 1;            % One-sided guard cells in velocity (columns)
        RankK = [];                   % OS-CFAR rank statistic
    end
    
    methods
        function obj = CfarDetector2D(varargin)
            % Construct CfarDetector2D instance with optional parameters.
            if nargin > 0
                for k = 1:2:nargin
                    propName = varargin{k};
                    propVal = varargin{k+1};
                    if isprop(obj, propName)
                        obj.(propName) = propVal;
                    else
                        error('CfarDetector2D:UnknownProperty', 'Property %s not recognized.', propName);
                    end
                end
            end
        end
        
        function [detectionMask, thresholdMatrix, noiseFloorMatrix] = detect(obj, rdmPower)
            % Execute 2D CFAR detection across a linear Range-Doppler power matrix.
            %
            % Inputs:
            %   rdmPower          - [N_range x N_doppler] non-negative power matrix
            %
            % Returns:
            %   detectionMask     - [N_range x N_doppler] boolean detection hit matrix
            %   thresholdMatrix   - [N_range x N_doppler] adaptive threshold values
            %   noiseFloorMatrix  - [N_range x N_doppler] estimated noise floor values
            
            [nR, nV] = size(rdmPower);
            
            trR = obj.TrainRange;
            gdR = obj.GuardRange;
            trV = obj.TrainVelocity;
            gdV = obj.GuardVelocity;
            
            % Outer and inner window dimensions
            outerDimR = 2 * (trR + gdR) + 1;
            outerDimV = 2 * (trV + gdV) + 1;
            guardDimR = 2 * gdR + 1;
            guardDimV = 2 * gdV + 1;
            
            nTotalTrain = (outerDimR * outerDimV) - (guardDimR * guardDimV);
            
            % Threshold scale factor alpha
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
            
            noiseFloorMatrix = zeros(nR, nV);
            thresholdMatrix = zeros(nR, nV);
            
            % Check if CA mode: use ultra-fast 2D integral image
            if strcmpi(obj.Mode, 'CA')
                % Pad matrix with zeros to handle 1-based indexing in integral lookups
                intImg = zeros(nR + 1, nV + 1);
                intImg(2:end, 2:end) = cumsum(cumsum(rdmPower, 1), 2);
                
                halfR = trR + gdR;
                halfV = trV + gdV;
                
                for r = 1:nR
                    r1_out = max(1, r - halfR);
                    r2_out = min(nR, r + halfR);
                    r1_in = max(1, r - gdR);
                    r2_in = min(nR, r + gdR);
                    
                    for v = 1:nV
                        v1_out = max(1, v - halfV);
                        v2_out = min(nV, v + halfV);
                        v1_in = max(1, v - gdV);
                        v2_in = min(nV, v + gdV);
                        
                        % Outer sum via integral image
                        sumOut = intImg(r2_out + 1, v2_out + 1) - ...
                                 intImg(r1_out, v2_out + 1) - ...
                                 intImg(r2_out + 1, v1_out) + ...
                                 intImg(r1_out, v1_out);
                                 
                        % Inner (guard + CUT) sum
                        sumIn = intImg(r2_in + 1, v2_in + 1) - ...
                                intImg(r1_in, v2_in + 1) - ...
                                intImg(r2_in + 1, v1_in) + ...
                                intImg(r1_in, v1_in);
                                
                        areaOut = (r2_out - r1_out + 1) * (v2_out - v1_out + 1);
                        areaIn = (r2_in - r1_in + 1) * (v2_in - v1_in + 1);
                        trainArea = areaOut - areaIn;
                        
                        if trainArea > 0
                            nf = (sumOut - sumIn) / trainArea;
                        else
                            nf = 1e-12;
                        end
                        
                        noiseFloorMatrix(r, v) = nf;
                        thresholdMatrix(r, v) = alpha * nf;
                    end
                end
            else
                % OS, GO, SO modes: extract training cells directly
                halfR = trR + gdR;
                halfV = trV + gdV;
                
                for r = 1:nR
                    r1_out = max(1, r - halfR);
                    r2_out = min(nR, r + halfR);
                    r1_in = max(1, r - gdR);
                    r2_in = min(nR, r + gdR);
                    
                    for v = 1:nV
                        v1_out = max(1, v - halfV);
                        v2_out = min(nV, v + halfV);
                        v1_in = max(1, v - gdV);
                        v2_in = min(nV, v + gdV);
                        
                        outerWindow = rdmPower(r1_out:r2_out, v1_out:v2_out);
                        % Create mask of guard region inside outer window
                        guardMask = false(size(outerWindow));
                        grd_r1 = r1_in - r1_out + 1;
                        grd_r2 = r2_in - r1_out + 1;
                        grd_v1 = v1_in - v1_out + 1;
                        grd_v2 = v2_in - v1_out + 1;
                        guardMask(grd_r1:grd_r2, grd_v1:grd_v2) = true;
                        
                        trainCells = outerWindow(~guardMask);
                        
                        if isempty(trainCells)
                            nf = 1e-12;
                        else
                            switch upper(obj.Mode)
                                case 'OS'
                                    sortedCells = sort(trainCells);
                                    kIdx = min(kRank, length(sortedCells));
                                    nf = sortedCells(kIdx);
                                case 'GO'
                                    halfIdx = floor(length(trainCells) / 2);
                                    nf = max(mean(trainCells(1:halfIdx)), mean(trainCells(halfIdx+1:end)));
                                case 'SO'
                                    halfIdx = floor(length(trainCells) / 2);
                                    nf = min(mean(trainCells(1:halfIdx)), mean(trainCells(halfIdx+1:end)));
                                otherwise
                                    nf = mean(trainCells);
                            end
                        end
                        
                        noiseFloorMatrix(r, v) = nf;
                        thresholdMatrix(r, v) = alpha * nf;
                    end
                end
            end
            
            detectionMask = (rdmPower > thresholdMatrix);
        end
    end
end
