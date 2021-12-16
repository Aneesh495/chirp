classdef IqCalibrator
    % IQCALIBRATOR Receiver baseband I/Q amplitude and phase imbalance compensation.
    %
    % Estimates and corrects DC offsets, quadrature phase error (deviation from 90 deg),
    % and gain mismatch between in-phase and quadrature ADC acquisition channels.
    
    properties
        DcOffsetI = 0.0;
        DcOffsetQ = 0.0;
        GainMismatch = 1.0;     % Ratio A_I / A_Q
        PhaseErrorRad = 0.0;    % Orthogonality error in radians
    end
    
    methods
        function obj = IqCalibrator()
            % Construct an uncalibrated default instance.
        end
        
        function obj = estimateBlind(obj, rawIqSignal)
            % Estimate I/Q imbalance parameters using Gram-Schmidt blind moments.
            %
            % Model:
            %   I_raw = I_ideal + dcI
            %   Q_raw = gain * (sin(phi)*I_ideal + cos(phi)*Q_ideal) + dcQ
            
            rawIq = rawIqSignal(:);
            I = real(rawIq);
            Q = imag(rawIq);
            
            % Step 1: DC offset removal
            obj.DcOffsetI = mean(I);
            obj.DcOffsetQ = mean(Q);
            
            I_zeroMean = I - obj.DcOffsetI;
            Q_zeroMean = Q - obj.DcOffsetQ;
            
            % Step 2: Power and correlation moments
            pI = mean(I_zeroMean.^2);
            pQ = mean(Q_zeroMean.^2);
            crossIQ = mean(I_zeroMean .* Q_zeroMean);
            
            % Step 3: Imbalance extraction
            % sin(phi) = crossIQ / (sqrt(pI) * sqrt(pQ))
            sinPhi = crossIQ / max(sqrt(pI * pQ), eps);
            sinPhi = max(-0.999, min(0.999, sinPhi));
            obj.PhaseErrorRad = asin(sinPhi);
            
            % Gain mismatch ratio
            obj.GainMismatch = sqrt(pI / max(pQ, eps));
        end
        
        function calibrated = correct(obj, rawIqSignal)
            % Apply inverse calibration matrix to raw received I/Q signal:
            %   [I_c; Q_c] = [1, 0; -tan(phi), 1/(gain * cos(phi))] * [I - dcI; Q - dcQ]
            
            I = real(rawIqSignal) - obj.DcOffsetI;
            Q = imag(rawIqSignal) - obj.DcOffsetQ;
            
            cosPhi = cos(obj.PhaseErrorRad);
            tanPhi = tan(obj.PhaseErrorRad);
            
            I_corr = I;
            Q_corr = (-tanPhi * I + (1.0 / (obj.GainMismatch * cosPhi)) * Q);
            
            calibrated = I_corr + 1j * Q_corr;
        end
    end
end
