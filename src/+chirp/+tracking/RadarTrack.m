classdef RadarTrack < handle
    % RADARTRACK State representation and tracking lifecycle metadata.
    %
    % Encapsulates 6D kinematic state [x; y; z; vx; vy; vz], error covariance P,
    % track status, confirmation counters, and trajectory history.
    
    properties
        TrackId                       % Unique integer track ID
        State = zeros(6, 1);          % 6x1 state vector [x; y; z; vx; vy; vz]
        Covariance = eye(6);          % 6x6 state estimation error covariance matrix P
        Status = 'tentative';         % Status: 'tentative', 'confirmed', 'coasted', 'deleted'
        Age = 1;                      % Total frames since track creation
        Hits = 1;                     % Total successful detection associations
        ConsecutiveMisses = 0;        % Number of consecutive frames without association
        HistoryPositions = zeros(3, 0); % Trajectory centroid history [3 x N_hist]
        LastTimestamp = 0.0;          % Last update timestamp (s)
    end
    
    methods
        function obj = RadarTrack(trackId, initialState, initialCov, timestamp)
            % Construct a RadarTrack instance.
            if nargin >= 1, obj.TrackId = trackId; end
            if nargin >= 2 && ~isempty(initialState), obj.State = initialState(:); end
            if nargin >= 3 && ~isempty(initialCov), obj.Covariance = initialCov; end
            if nargin >= 4, obj.LastTimestamp = timestamp; end
            
            if ~isempty(obj.State)
                obj.HistoryPositions = obj.State(1:3);
            end
        end
        
        function updateHistory(obj)
            % Append current position to track trajectory history.
            obj.HistoryPositions(:, end + 1) = obj.State(1:3);
            if size(obj.HistoryPositions, 2) > 100
                obj.HistoryPositions(:, 1) = []; % Retain up to 100 recent steps
            end
        end
        
        function pos = getPosition(obj)
            pos = obj.State(1:3);
        end
        
        function vel = getVelocity(obj)
            vel = obj.State(4:6);
        end
        
        function speed = getSpeed(obj)
            speed = sqrt(sum(obj.State(4:6).^2));
        end
    end
end
