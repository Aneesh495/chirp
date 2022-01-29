classdef TrackLifecycleManager
    % TRACKLIFECYCLEMANAGER Track quality control, confirmation, and duplicate merging.
    %
    % Manages M-out-of-N track confirmation policies, coasting inflation, and
    % merges duplicate tracks occupying overlapping spatial ellipsoids.
    
    properties
        M_Confirm = 3;                % Hits required within N window
        N_Window = 5;                 % Sliding window length
        MaxCoastFrames = 3;           % Maximum consecutive missed frames before deletion
        MergeDistanceThreshold = 1.5; % Mahalanobis distance threshold for merging duplicate tracks
    end
    
    methods
        function obj = TrackLifecycleManager(mConfirm, nWindow, maxCoast)
            % Construct a TrackLifecycleManager instance.
            if nargin >= 1 && ~isempty(mConfirm), obj.M_Confirm = mConfirm; end
            if nargin >= 2 && ~isempty(nWindow), obj.N_Window = nWindow; end
            if nargin >= 3 && ~isempty(maxCoast), obj.MaxCoastFrames = maxCoast; end
        end
        
        function cleanedTracks = pruneAndMerge(obj, tracks)
            % Prune dead tracks and merge redundant duplicate tracks.
            %
            % Input:
            %   tracks        - Array of chirp.tracking.RadarTrack objects
            % Output:
            %   cleanedTracks - Filtered array of distinct active tracks
            
            if isempty(tracks)
                cleanedTracks = tracks;
                return;
            end
            
            % Step 1: Remove dead tracks
            aliveMask = true(1, length(tracks));
            for i = 1:length(tracks)
                if tracks(i).ConsecutiveMisses > obj.MaxCoastFrames
                    aliveMask(i) = false;
                end
            end
            tracks = tracks(aliveMask);
            
            % Step 2: Merge duplicate tracks
            n = length(tracks);
            if n <= 1
                cleanedTracks = tracks;
                return;
            end
            
            toDelete = false(1, n);
            for i = 1:n
                if toDelete(i), continue; end
                p1 = tracks(i).State(1:3);
                P1 = tracks(i).Covariance(1:3, 1:3);
                
                for j = (i + 1):n
                    if toDelete(j), continue; end
                    p2 = tracks(j).State(1:3);
                    P2 = tracks(j).Covariance(1:3, 1:3);
                    
                    diffP = p1 - p2;
                    combinedCov = P1 + P2;
                    distM = sqrt(diffP' * (combinedCov \ diffP));
                    
                    if distM < obj.MergeDistanceThreshold
                        % Merge j into i: retain track with higher hit count
                        if tracks(j).Hits > tracks(i).Hits
                            tracks(i).State = tracks(j).State;
                            tracks(i).Covariance = tracks(j).Covariance;
                            tracks(i).Hits = tracks(j).Hits;
                        end
                        toDelete(j) = true;
                    end
                end
            end
            
            cleanedTracks = tracks(~toDelete);
        end
    end
end
