classdef ExtendedKalmanFilter < handle
    % EXTENDEDKALMANFILTER Non-linear radar kinematic state estimator.
    %
    % Estimates 3D Cartesian position and velocity [x; y; z; vx; vy; vz] from
    % spherical radar measurements: range r, azimuth az, elevation el, and Doppler vr.
    
    properties
        ProcessNoiseSigma = 0.8;      % Process acceleration noise standard deviation (m/s^2)
        MeasurementNoiseR             % 4x4 measurement noise covariance R [r, az(rad), el(rad), vr]
    end
    
    methods
        function obj = ExtendedKalmanFilter(measNoiseR, procNoiseSigma)
            % Construct an ExtendedKalmanFilter instance.
            if nargin >= 1 && ~isempty(measNoiseR)
                obj.MeasurementNoiseR = measNoiseR;
            else
                % Default radar measurement variances:
                % sigma_r = 0.25 m, sigma_az = 0.5 deg, sigma_el = 1.0 deg, sigma_vr = 0.2 m/s
                varR = 0.25^2;
                varAz = (0.5 * pi / 180.0)^2;
                varEl = (1.0 * pi / 180.0)^2;
                varVr = 0.2^2;
                obj.MeasurementNoiseR = diag([varR, varAz, varEl, varVr]);
            end
            if nargin >= 2 && ~isempty(procNoiseSigma)
                obj.ProcessNoiseSigma = procNoiseSigma;
            end
        end
        
        function [predState, predCov] = predict(obj, state, cov, dt)
            % Propagate state and error covariance forward by dt using Constant Velocity model.
            [F, Q] = chirp.tracking.MotionModel.constantVelocity(dt, obj.ProcessNoiseSigma);
            predState = F * state;
            predCov = F * cov * F' + Q;
        end
        
        function [zPred, H] = measurementModel(obj, state)
            % Non-linear observation function h(x) and analytical Jacobian H.
            %
            % State: [x; y; z; vx; vy; vz]
            % Measurement: [r; az(rad); el(rad); vr]
            
            x = state(1); y = state(2); z = state(3);
            vx = state(4); vy = state(5); vz = state(6);
            
            rSq = x^2 + y^2 + z^2;
            r = sqrt(max(rSq, 1e-6));
            rxySq = x^2 + y^2;
            rxy = sqrt(max(rxySq, 1e-6));
            
            az = atan2(x, y);
            el = atan2(z, rxy);
            vr = (x * vx + y * vy + z * vz) / r;
            
            zPred = [r; az; el; vr];
            
            % Compute 4x6 Jacobian matrix H
            H = zeros(4, 6);
            
            % dr / dx, y, z
            H(1, 1) = x / r;
            H(1, 2) = y / r;
            H(1, 3) = z / r;
            
            % daz / dx, y, z (az = atan2(x, y))
            H(2, 1) = y / rxySq;
            H(2, 2) = -x / rxySq;
            
            % del / dx, y, z (el = atan2(z, rxy))
            H(3, 1) = -x * z / (rSq * rxy);
            H(3, 2) = -y * z / (rSq * rxy);
            H(3, 3) = rxy / rSq;
            
            % dvr / dx, y, z, vx, vy, vz
            H(4, 1) = vx / r - x * vr / rSq;
            H(4, 2) = vy / r - y * vr / rSq;
            H(4, 3) = vz / r - z * vr / rSq;
            H(4, 4) = x / r;
            H(4, 5) = y / r;
            H(4, 6) = z / r;
        end
        
        function [updState, updCov, innov, S] = update(obj, predState, predCov, measurement)
            % Update predicted state with measurement z = [r; az(rad); el(rad); vr].
            
            zMeas = measurement(:);
            [zPred, H] = obj.measurementModel(predState);
            
            % Innovation vector
            innov = zMeas - zPred;
            % Wrap angular innovation in radians to [-pi, pi)
            innov(2) = chirp.utils.FastMath.wrap_to_pi(innov(2));
            innov(3) = chirp.utils.FastMath.wrap_to_pi(innov(3));
            
            % Innovation covariance S = H * P * H' + R
            R = obj.MeasurementNoiseR;
            S = H * predCov * H' + R;
            
            % Kalman gain: K = P * H' * S^(-1)
            K = (predCov * H') / S;
            
            % State update
            updState = predState + K * innov;
            
            % Joseph-form covariance update for numerical symmetry:
            % P = (I - K*H) * P * (I - K*H)' + K * R * K'
            I6 = eye(6);
            IKH = I6 - K * H;
            updCov = IKH * predCov * IKH' + K * R * K';
        end
    end
end
