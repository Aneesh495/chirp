classdef GnnTracker < handle
    % GNNTRACKER Global Nearest Neighbor multi-target kinematic radar tracker.
    %
    % Performs statistical gating, global bipartite Hungarian assignment,
    % EKF track filtering, and track lifecycle management.
    
    properties
        Tracks = [];                  % Array of active RadarTrack objects
        NextTrackId = 1;              % Counter for track ID allocation
        Ekf                           % ExtendedKalmanFilter instance
        GateThreshold = 16.0;         % Mahalanobis distance chi-square gating threshold (99.7%)
        ConfirmationThreshold = 3;    % Minimum hits to confirm track
        MaxConsecutiveMisses = 3;     % Misses before track deletion
    end
    
    methods
        function obj = GnnTracker(ekfInstance, gateThreshold)
            % Construct a GnnTracker instance.
            if nargin >= 1 && ~isempty(ekfInstance)
                obj.Ekf = ekfInstance;
            else
                obj.Ekf = chirp.tracking.ExtendedKalmanFilter();
            end
            if nargin >= 2 && ~isempty(gateThreshold)
                obj.GateThreshold = gateThreshold;
            end
        end
        
        function activeTracks = step(obj, detections, timestamp, dt)
            % Execute one tracking cycle for incoming radar detections [4 x N_det].
            %
            % Inputs:
            %   detections - [4 x N_det] matrix: [range; az(rad); el(rad); vr]
            %   timestamp  - Current frame epoch timestamp (s)
            %   dt         - Time step elapsed since last frame (s)
            %
            % Returns:
            %   activeTracks - Array of active confirmed / tentative RadarTrack instances
            
            nTracks = length(obj.Tracks);
            nDet = size(detections, 2);
            
            % Step 1: Predict all tracks forward
            for iTrk = 1:nTracks
                trk = obj.Tracks(iTrk);
                [pState, pCov] = obj.Ekf.predict(trk.State, trk.Covariance, dt);
                trk.State = pState;
                trk.Covariance = pCov;
                trk.Age = trk.Age + 1;
            end
            
            % Step 2: Compute Mahalanobis cost matrix
            costMatrix = Inf(nTracks, nDet);
            
            for iTrk = 1:nTracks
                trk = obj.Tracks(iTrk);
                [zPred, H] = obj.Ekf.measurementModel(trk.State);
                S = H * trk.Covariance * H' + obj.Ekf.MeasurementNoiseR;
                invS = inv(S);
                
                for jDet = 1:nDet
                    zMeas = detections(:, jDet);
                    innov = zMeas - zPred;
                    innov(2) = chirp.utils.FastMath.wrap_to_pi(innov(2));
                    innov(3) = chirp.utils.FastMath.wrap_to_pi(innov(3));
                    
                    mahalanobisSq = real(innov' * invS * innov);
                    if mahalanobisSq <= obj.GateThreshold
                        costMatrix(iTrk, jDet) = mahalanobisSq;
                    end
                end
            end
            
            % Step 3: Optimal assignment
            [matches, unassignedTracks, unassignedDets] = ...
                chirp.tracking.HungarianMatcher.match(costMatrix, obj.GateThreshold);
                
            % Step 4: Update matched tracks
            if ~isempty(matches)
                for k = 1:size(matches, 1)
                    iTrk = matches(k, 1);
                    jDet = matches(k, 2);
                    trk = obj.Tracks(iTrk);
                    
                    zMeas = detections(:, jDet);
                    [uState, uCov, ~, ~] = obj.Ekf.update(trk.State, trk.Covariance, zMeas);
                    
                    trk.State = uState;
                    trk.Covariance = uCov;
                    trk.Hits = trk.Hits + 1;
                    trk.ConsecutiveMisses = 0;
                    trk.LastTimestamp = timestamp;
                    trk.updateHistory();
                    
                    if trk.Hits >= obj.ConfirmationThreshold
                        trk.Status = 'confirmed';
                    end
                end
            end
            
            % Step 5: Coast unassigned tracks
            for k = 1:length(unassignedTracks)
                iTrk = unassignedTracks(k);
                trk = obj.Tracks(iTrk);
                trk.ConsecutiveMisses = trk.ConsecutiveMisses + 1;
                trk.Status = 'coasted';
            end
            
            % Step 6: Initialize new tentative tracks from unassigned detections
            for k = 1:length(unassignedDets)
                jDet = unassignedDets(k);
                z = detections(:, jDet);
                r = z(1); az = z(2); el = z(3); vr = z(4);
                
                % Convert spherical observation to initial Cartesian state
                x = r * cos(el) * sin(az);
                y = r * cos(el) * cos(az);
                zCoord = r * sin(el);
                
                % Velocity along line-of-sight
                vx = vr * sin(az) * cos(el);
                vy = vr * cos(az) * cos(el);
                vz = vr * sin(el);
                
                initState = [x; y; zCoord; vx; vy; vz];
                initCov = diag([1.0, 1.0, 1.0, 4.0, 4.0, 4.0]);
                
                newTrk = chirp.tracking.RadarTrack(obj.NextTrackId, initState, initCov, timestamp);
                obj.NextTrackId = obj.NextTrackId + 1;
                if isempty(obj.Tracks)
                    obj.Tracks = newTrk;
                else
                    obj.Tracks(end + 1) = newTrk;
                end
            end
            
            % Step 7: Delete tracks exceeding consecutive miss threshold
            if ~isempty(obj.Tracks)
                survivingMask = true(1, length(obj.Tracks));
                for iTrk = 1:length(obj.Tracks)
                    if obj.Tracks(iTrk).ConsecutiveMisses > obj.MaxConsecutiveMisses
                        survivingMask(iTrk) = false;
                    end
                end
                obj.Tracks = obj.Tracks(survivingMask);
            end
            
            activeTracks = obj.Tracks;
        end
    end
end
