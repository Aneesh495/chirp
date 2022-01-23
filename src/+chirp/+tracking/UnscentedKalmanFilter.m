classdef UnscentedKalmanFilter < handle
    % UNSCENTEDKALMANFILTER Non-linear radar state estimator via Scaled Unscented Transform.
    %
    % Propagates deterministic sigma points through the non-linear radar observation
    % function to capture posterior mean and covariance without truncation errors.
    
    properties
        Alpha = 1e-3;                 % Primary spread parameter (0 < alpha <= 1)
        Beta = 2.0;                   % Prior distribution parameter (beta = 2 optimal for Gaussian)
        Kappa = 0.0;                  % Secondary scaling parameter
        ProcessNoiseSigma = 0.8;      % Process noise intensity
        MeasurementNoiseR             % 4x4 measurement noise covariance R
    end
    
    methods
        function obj = UnscentedKalmanFilter(measNoiseR, procNoiseSigma)
            % Construct an UnscentedKalmanFilter instance.
            if nargin >= 1 && ~isempty(measNoiseR)
                obj.MeasurementNoiseR = measNoiseR;
            else
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
            % Constant velocity prediction step.
            [F, Q] = chirp.tracking.MotionModel.constantVelocity(dt, obj.ProcessNoiseSigma);
            predState = F * state;
            predCov = F * cov * F' + Q;
        end
        
        function [sigmaPoints, Wm, Wc] = generateSigmaPoints(obj, state, cov)
            % Compute 2L + 1 sigma points and associated mean/covariance weights.
            L = length(state);
            lambdaVal = (obj.Alpha^2) * (L + obj.Kappa) - L;
            gammaVal = sqrt(L + lambdaVal);
            
            % Matrix square root via Cholesky
            % Ensure positive definite
            covSym = 0.5 * (cov + cov');
            [S_chol, p] = chol(covSym, 'lower');
            if p ~= 0
                % Regularize if ill-conditioned
                S_chol = chol(covSym + 1e-6 * eye(L), 'lower');
            end
            
            sigmaPoints = zeros(L, 2 * L + 1);
            sigmaPoints(:, 1) = state;
            
            for i = 1:L
                col = S_chol(:, i);
                sigmaPoints(:, 1 + i) = state + gammaVal * col;
                sigmaPoints(:, 1 + L + i) = state - gammaVal * col;
            end
            
            % Compute weights
            Wm = zeros(1, 2 * L + 1);
            Wc = zeros(1, 2 * L + 1);
            
            Wm(1) = lambdaVal / (L + lambdaVal);
            Wc(1) = Wm(1) + (1.0 - obj.Alpha^2 + obj.Beta);
            
            w_i = 1.0 / (2.0 * (L + lambdaVal));
            Wm(2:end) = w_i;
            Wc(2:end) = w_i;
        end
        
        function [updState, updCov, innov, S] = update(obj, predState, predCov, measurement)
            % Non-linear observation update using Unscented Transform.
            
            zMeas = measurement(:);
            L = length(predState);
            [sigmaPts, Wm, Wc] = obj.generateSigmaPoints(predState, predCov);
            
            numPts = 2 * L + 1;
            Z_sigma = zeros(4, numPts);
            
            % Transform each sigma point through observation function h(x)
            for i = 1:numPts
                pt = sigmaPts(:, i);
                x = pt(1); y = pt(2); z = pt(3);
                vx = pt(4); vy = pt(5); vz = pt(6);
                
                r = sqrt(max(x^2 + y^2 + z^2, 1e-6));
                rxy = sqrt(max(x^2 + y^2, 1e-6));
                az = atan2(x, y);
                el = atan2(z, rxy);
                vr = (x * vx + y * vy + z * vz) / r;
                
                Z_sigma(:, i) = [r; az; el; vr];
            end
            
            % Predicted measurement mean z_hat = sum(Wm_i * Z_i)
            zPred = Z_sigma * Wm.';
            
            % Measurement covariance S and cross-covariance Pxz
            S = obj.MeasurementNoiseR;
            Pxz = zeros(L, 4);
            
            for i = 1:numPts
                diffZ = Z_sigma(:, i) - zPred;
                diffZ(2) = chirp.utils.FastMath.wrap_to_pi(diffZ(2));
                diffZ(3) = chirp.utils.FastMath.wrap_to_pi(diffZ(3));
                
                diffX = sigmaPts(:, i) - predState;
                
                S = S + Wc(i) * (diffZ * diffZ');
                Pxz = Pxz + Wc(i) * (diffX * diffZ');
            end
            
            % Kalman gain K = Pxz * S^(-1)
            K = Pxz / S;
            
            % Innovation
            innov = zMeas - zPred;
            innov(2) = chirp.utils.FastMath.wrap_to_pi(innov(2));
            innov(3) = chirp.utils.FastMath.wrap_to_pi(innov(3));
            
            updState = predState + K * innov;
            updCov = predCov - K * S * K';
            updCov = 0.5 * (updCov + updCov'); % Ensure symmetry
        end
    end
end
