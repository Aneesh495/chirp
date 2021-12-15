classdef RangeDopplerMap
    % RANGEDOPPLERMAP Container and feature extractor for 2D Range-Doppler representations.
    %
    % Manages complex and power spectral data, non-coherent channel integration,
    % dynamic range compression, thresholding, and peak detection.
    
    properties
        Data                          % Complex [N_range x N_doppler x N_virt] or power [N_range x N_doppler]
        RangeAxis                     % 1 x N_range vector of range coordinates (m)
        VelocityAxis                  % 1 x N_doppler vector of radial velocity coordinates (m/s)
    end
    
    methods
        function obj = RangeDopplerMap(data, rangeAxis, velocityAxis)
            % Construct a RangeDopplerMap container.
            if nargin >= 3
                obj.Data = data;
                obj.RangeAxis = rangeAxis(:).';
                obj.VelocityAxis = velocityAxis(:).';
            end
        end
        
        function pwrDb = nonCoherentIntegration(obj)
            % Compute Non-Coherent Integration (NCI) power matrix in decibels:
            %   P(r, v) = 10 * log10( sum_k |Data(r, v, k)|^2 / N_virt )
            
            if ndims(obj.Data) == 3
                pwrLin = mean(abs(obj.Data).^2, 3);
            else
                pwrLin = abs(obj.Data).^2;
            end
            pwrDb = 10.0 * log10(max(pwrLin, eps));
        end
        
        function pwrDb = coherentIntegration(obj, steerVec)
            % Compute Coherent Integration power matrix focused along a steering vector.
            %
            % Input:
            %   steerVec - [N_virt x 1] complex unit norm steering vector
            %
            % Output:
            %   pwrDb    - [N_range x N_doppler] power in decibels
            
            if ndims(obj.Data) ~= 3
                error('RangeDopplerMap:MultiChannelRequired', 'Coherent beamforming requires multi-channel 3D cube.');
            end
            
            [nR, nV, nVirt] = size(obj.Data);
            sv = steerVec(:);
            if length(sv) ~= nVirt
                error('RangeDopplerMap:SteeringVectorMismatch', 'Steering vector length must match virtual channel count.');
            end
            
            % Matrix multiplication across channels
            reshaped = reshape(obj.Data, [nR * nV, nVirt]);
            beamformed = reshaped * conj(sv);
            beamMatrix = reshape(beamformed, [nR, nV]);
            
            pwrDb = 10.0 * log10(max(abs(beamMatrix).^2, eps));
        end
        
        function [pks, rangeCoords, velCoords] = findPeaks(obj, thresholdDb, minDistanceBins)
            % Find local 2D spectral peaks exceeding thresholdDb with minimum bin separation.
            
            pwr = obj.nonCoherentIntegration();
            if nargin < 2 || isempty(thresholdDb)
                thresholdDb = max(pwr(:)) - 25.0;
            end
            if nargin < 3 || isempty(minDistanceBins)
                minDistanceBins = 3;
            end
            
            [nR, nV] = size(pwr);
            
            % Local maximum test with minimum distance dilation
            peaksMask = false(nR, nV);
            d = minDistanceBins;
            
            for r = (d + 1):(nR - d)
                for v = (d + 1):(nV - d)
                    val = pwr(r, v);
                    if val >= thresholdDb
                        neighborhood = pwr((r - d):(r + d), (v - d):(v + d));
                        % Exclude center cell to enforce strict local peak
                        neighborhood(d + 1, d + 1) = -Inf;
                        if val > max(neighborhood(:))
                            peaksMask(r, v) = true;
                        end
                    end
                end
            end
            
            [rIdx, vIdx] = find(peaksMask);
            if isempty(rIdx)
                pks = [];
                rangeCoords = [];
                velCoords = [];
                return;
            end
            
            pks = zeros(length(rIdx), 1);
            rangeCoords = zeros(length(rIdx), 1);
            velCoords = zeros(length(rIdx), 1);
            
            for k = 1:length(rIdx)
                pks(k) = pwr(rIdx(k), vIdx(k));
                rangeCoords(k) = obj.RangeAxis(rIdx(k));
                velCoords(k) = obj.VelocityAxis(vIdx(k));
            end
            
            % Sort detected peaks in descending order of power
            [pks, sortIdx] = sort(pks, 'descend');
            rangeCoords = rangeCoords(sortIdx);
            velCoords = velCoords(sortIdx);
        end
    end
end
