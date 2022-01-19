classdef MotionModel
    % MOTIONMODEL Kinematic state transition and process noise covariance matrices.
    %
    % Formulates discrete state transition matrices F and process noise Q for
    % Constant Velocity (CV) and Constant Acceleration (CA) Newtonian motion models.
    
    methods (Static)
        function [F, Q] = constantVelocity(dt, processNoiseSigma)
            % 6-State Constant Velocity (CV) model: [x; y; z; vx; vy; vz]
            %
            % State transition:
            %   p(k+1) = p(k) + v(k) * dt
            %   v(k+1) = v(k)
            %
            % Continuous white noise acceleration model Q_k.
            
            if nargin < 2 || isempty(processNoiseSigma)
                processNoiseSigma = 0.5; % m/s^2 acceleration noise intensity
            end
            
            I3 = eye(3);
            O3 = zeros(3);
            
            F = [I3, dt * I3; ...
                 O3, I3];
                 
            % Discrete continuous white noise acceleration covariance:
            % Q_11 = (dt^3 / 3) * sigma^2 * I3
            % Q_12 = Q_21 = (dt^2 / 2) * sigma^2 * I3
            % Q_22 = dt * sigma^2 * I3
            q11 = (dt^3 / 3.0) * (processNoiseSigma^2) * I3;
            q12 = (dt^2 / 2.0) * (processNoiseSigma^2) * I3;
            q22 = dt * (processNoiseSigma^2) * I3;
            
            Q = [q11, q12; ...
                 q12, q22];
        end
        
        function [F, Q] = constantAcceleration(dt, jerkNoiseSigma)
            % 9-State Constant Acceleration (CA) model: [p; v; a]
            if nargin < 2 || isempty(jerkNoiseSigma)
                jerkNoiseSigma = 0.2; % m/s^3 jerk noise intensity
            end
            
            I3 = eye(3);
            O3 = zeros(3);
            
            F = [I3, dt * I3, 0.5 * (dt^2) * I3; ...
                 O3, I3,      dt * I3; ...
                 O3, O3,      I3];
                 
            q11 = (dt^5 / 20.0) * (jerkNoiseSigma^2) * I3;
            q12 = (dt^4 / 8.0) * (jerkNoiseSigma^2) * I3;
            q13 = (dt^3 / 6.0) * (jerkNoiseSigma^2) * I3;
            q22 = (dt^3 / 3.0) * (jerkNoiseSigma^2) * I3;
            q23 = (dt^2 / 2.0) * (jerkNoiseSigma^2) * I3;
            q33 = dt * (jerkNoiseSigma^2) * I3;
            
            Q = [q11, q12, q13; ...
                 q12, q22, q23; ...
                 q13, q23, q33];
        end
    end
end
