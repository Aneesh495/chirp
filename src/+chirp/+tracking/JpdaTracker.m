classdef JpdaTracker < handle
    % JPDATRACKER Joint Probabilistic Data Association multi-target radar tracker.
    %
    % Performs soft probabilistic association across overlapping validation gates
    % to maintain robust track continuity in high-density clutter environments.
    
    properties
        Tracks = [];
        NextTrackId = 1;
        Ekf
        ClutterDensity = 1e-4;        % Spatial clutter rate lambda_c
        GateThreshold = 16.0;         % Validation gate threshold
    end
    
    methods
        function obj = JpdaTracker(ekfInstance, clutterDensity)
            % Construct a JpdaTracker instance.
            if nargin >= 1 && ~isempty(ekfInstance)
                obj.Ekf = ekfInstance;
            else
                obj.Ekf = chirp.tracking.ExtendedKalmanFilter();
            end
            if nargin >= 2 && ~isempty(clutterDensity)
                obj.ClutterDensity = clutterDensity;
            end
        end
        
        function activeTracks = step(obj, detections, timestamp, dt)
            % Execute one JPDA tracking update cycle.
            nTracks = length(obj.Tracks);
            nDet = size(detections, 2);
            
            % Step 1: Predict tracks
            for iTrk = 1:nTracks
                trk = obj.Tracks(iTrk);
                [pState, pCov] = obj.Ekf.predict(trk.State, trk.Covariance, dt);
                trk.State = pState;
                trk.Covariance = pCov;
                trk.Age = trk.Age + 1;
            end
            
            if nTracks == 0
                % Initialize tracks from all detections
                for j = 1:nDet
                    z = detections(:, j);
                    r = z(1); az = z(2); el = z(3); vr = z(4);
                    x = r * cos(el) * sin(az);
                    y = r * cos(el) * cos(az);
                    zCoord = r * sin(el);
                    vx = vr * sin(az) * cos(el);
                    vy = vr * cos(az) * cos(el);
                    vz = vr * sin(el);
                    
                    initState = [x; y; zCoord; vx; vy; vz];
                    initCov = diag([1.0, 1.0, 1.0, 4.0, 4.0, 4.0]);
                    newTrk = chirp.tracking.RadarTrack(...
                        obj.NextTrackId, initState, initCov, timestamp);
                    obj.NextTrackId = obj.NextTrackId + 1;
                    if isempty(obj.Tracks)
                        obj.Tracks = newTrk;
                    else
                        obj.Tracks(end + 1) = newTrk;
                    end
                end
                activeTracks = obj.Tracks;
                return;
            end
            
            % Step 2: Compute marginal association probabilities for each track
            for iTrk = 1:nTracks
                trk = obj.Tracks(iTrk);
                [zPred, H] = obj.Ekf.measurementModel(trk.State);
                S = H * trk.Covariance * H' + obj.Ekf.MeasurementNoiseR;
                invS = inv(S);
                detS = max(det(S), 1e-12);
                K = (trk.Covariance * H') / S;
                
                innovations = zeros(4, nDet);
                likelihoods = zeros(1, nDet);
                validMask = false(1, nDet);
                
                for jDet = 1:nDet
                    zMeas = detections(:, jDet);
                    y = zMeas - zPred;
                    y(2) = chirp.utils.FastMath.wrap_to_pi(y(2));
                    y(3) = chirp.utils.FastMath.wrap_to_pi(y(3));
                    innovations(:, jDet) = y;
                    
                    d2 = real(y' * invS * y);
                    if d2 <= obj.GateThreshold
                        validMask(jDet) = true;
                        normFactor = 1.0 / sqrt(((2.0 * pi)^4) * detS);
                        likelihoods(jDet) = normFactor * exp(-0.5 * d2);
                    end
                end
                
                sumL = sum(likelihoods(validMask));
                denom = obj.ClutterDensity + sumL;
                
                % Beta_0: probability of missed detection / clutter
                beta0 = obj.ClutterDensity / max(denom, 1e-12);
                
                % Combined innovation
                combInnov = zeros(4, 1);
                if sumL > 0
                    for jDet = find(validMask)
                        beta_j = likelihoods(jDet) / denom;
                        combInnov = combInnov + beta_j * innovations(:, jDet);
                    end
                    
                    % Update track state with combined innovation
                    trk.State = trk.State + K * combInnov;
                    trk.Hits = trk.Hits + 1;
                    trk.ConsecutiveMisses = 0;
                    trk.Status = 'confirmed';
                else
                    trk.ConsecutiveMisses = trk.ConsecutiveMisses + 1;
                end
                
                trk.LastTimestamp = timestamp;
                trk.updateHistory();
            end
            
            % Delete tracks exceeding miss threshold
            if ~isempty(obj.Tracks)
                survivingMask = true(1, length(obj.Tracks));
                for iTrk = 1:length(obj.Tracks)
                    if obj.Tracks(iTrk).ConsecutiveMisses > 4
                        survivingMask(iTrk) = false;
                    end
                end
                obj.Tracks = obj.Tracks(survivingMask);
            end
            
            activeTracks = obj.Tracks;
        end
    end
end
