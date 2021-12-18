classdef MutualCouplingCalibrator
    % MUTUALCOUPLINGCALIBRATOR Antenna array mutual electromagnetic coupling calibration.
    %
    % Models and decouples the mutual impedance coupling matrix across physical
    % and virtual receiver antenna elements:
    %   s_meas = C * s_ideal + n  =>  s_cal = C^(-1) * s_meas
    
    properties
        CouplingMatrix                % M x M complex coupling matrix C
        DecouplingMatrix              % M x M inverted decoupling matrix C^(-1)
    end
    
    methods
        function obj = MutualCouplingCalibrator(numElements, couplingCoeff)
            % Construct calibrator with Toeplitz symmetric coupling model.
            %
            % Inputs:
            %   numElements   - Integer number of antenna channels M
            %   couplingCoeff - Vector of off-diagonal coupling strengths [c1, c2, ...]
            
            if nargin >= 1
                if nargin < 2 || isempty(couplingCoeff)
                    couplingCoeff = [0.15 * exp(1j * pi / 4), 0.03 * exp(1j * pi / 6)];
                end
                
                cRow = zeros(1, numElements);
                cRow(1) = 1.0;
                kMax = min(length(couplingCoeff), numElements - 1);
                cRow(2:(kMax + 1)) = couplingCoeff(1:kMax);
                
                % Form symmetric Toeplitz coupling matrix
                C = toeplitz(cRow);
                obj.CouplingMatrix = C;
                obj.DecouplingMatrix = inv(C);
            end
        end
        
        function calibratedCube = decoupleCube(obj, rawCube)
            % Apply decoupling matrix across virtual antenna dimension (dim 3).
            %
            % Input:
            %   rawCube        - [N_range x N_doppler x N_virt]
            % Output:
            %   calibratedCube - Decoupled data cube with restored phase/amplitude
            
            [nR, nV, nVirt] = size(rawCube);
            if size(obj.DecouplingMatrix, 1) ~= nVirt
                error('MutualCouplingCalibrator:DimensionMismatch', ...
                    'Decoupling matrix dimension %d does not match virtual channels %d.', ...
                    size(obj.DecouplingMatrix, 1), nVirt);
            end
            
            reshaped = reshape(rawCube, [nR * nV, nVirt]);
            % s_cal = s_meas * (C^(-1)).'
            decoupled = reshaped * (obj.DecouplingMatrix.');
            calibratedCube = reshape(decoupled, [nR, nV, nVirt]);
        end
    end
end
