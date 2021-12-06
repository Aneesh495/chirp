classdef MimoArray
    % MIMOARRAY Virtual antenna array synthesis and steering vector generator.
    %
    % Computes spatial virtual element coordinates from physical TX and RX
    % element distributions, and synthesizes 3D spatial steering vectors.
    
    properties
        TxPositions = zeros(3, 1);    % 3 x N_tx physical transmit antenna coordinates (m)
        RxPositions = zeros(3, 4);    % 3 x N_rx physical receive antenna coordinates (m)
    end
    
    properties (Dependent)
        NumTx                         % Number of physical transmit antennas
        NumRx                         % Number of physical receive antennas
        NumVirtual                    % Number of virtual receive elements (N_tx * N_rx)
        VirtualPositions              % 3 x NumVirtual synthesized virtual antenna coordinates (m)
    end
    
    methods
        function obj = MimoArray(txPos, rxPos)
            % Construct a MIMO array from 3xN_tx and 3xN_rx antenna positions.
            if nargin >= 2
                if size(txPos, 1) ~= 3 || size(rxPos, 1) ~= 3
                    error('MimoArray:InvalidDimension', 'Antenna position arrays must have 3 rows [x; y; z].');
                end
                obj.TxPositions = txPos;
                obj.RxPositions = rxPos;
            end
        end
        
        function n = get.NumTx(obj)
            n = size(obj.TxPositions, 2);
        end
        
        function n = get.NumRx(obj)
            n = size(obj.RxPositions, 2);
        end
        
        function n = get.NumVirtual(obj)
            n = obj.NumTx * obj.NumRx;
        end
        
        function vPos = get.VirtualPositions(obj)
            % Virtual element position: r_v = r_tx + r_rx
            nTx = obj.NumTx;
            nRx = obj.NumRx;
            vPos = zeros(3, nTx * nRx);
            
            idx = 1;
            for iTx = 1:nTx
                tx = obj.TxPositions(:, iTx);
                for iRx = 1:nRx
                    rx = obj.RxPositions(:, iRx);
                    vPos(:, idx) = tx + rx;
                    idx = idx + 1;
                end
            end
        end
        
        function sv = steeringVector(obj, azDeg, elDeg, wavelengthMeters)
            % Compute spatial steering vector for incident plane wave from (azDeg, elDeg).
            %
            % Direction unit vector pointing from target towards radar (along boresight +y):
            %   kx = sin(az) * cos(el)
            %   ky = cos(az) * cos(el)
            %   kz = sin(el)
            %
            % Phase shift at virtual element position r_v:
            %   phi = (2*pi / lambda) * (r_v . k)
            
            azRad = azDeg * (pi / 180.0);
            elRad = elDeg * (pi / 180.0);
            
            kVec = [sin(azRad) * cos(elRad); ...
                    cos(azRad) * cos(elRad); ...
                    sin(elRad)];
            
            vPos = obj.VirtualPositions;
            % 1 x NumVirtual dot product
            spatialPhase = (2.0 * pi / wavelengthMeters) * (kVec.' * vPos);
            
            sv = exp(1j * spatialPhase(:));
        end
    end
    
    methods (Static)
        function arrayObj = automotive3Tx4Rx(wavelengthMeters)
            % Standard automotive 3 TX x 4 RX Time-Division Multiplexed (TDM) MIMO array.
            % Yields 12-element uniform linear virtual array spaced at lambda / 2 along X axis.
            
            d = wavelengthMeters / 2.0;
            
            % 4 RX antennas spaced by d: [0, d, 2d, 3d]
            rxPos = [ (0:3) * d; zeros(2, 4) ];
            
            % 3 TX antennas spaced by 4d: [0, 4d, 8d]
            txPos = [ [0, 4, 8] * d; zeros(2, 3) ];
            
            arrayObj = chirp.sim.MimoArray(txPos, rxPos);
        end
        
        function arrayObj = imaging12Tx16Rx(wavelengthMeters)
            % 4D Imaging Radar 12 TX x 16 RX MIMO array with elevation and azimuth baseline.
            % Yields 192 virtual channels.
            
            d = wavelengthMeters / 2.0;
            
            % 16 RX antennas in 2 rows of 8: spaced by d in X, row separation by d in Z
            [rxX, rxZ] = meshgrid((0:7) * d, [0, 1] * d);
            rxPos = [rxX(:).'; zeros(1, 16); rxZ(:).'];
            
            % 12 TX antennas in 3 rows of 4: spaced by 8d in X, 2d in Z
            [txX, txZ] = meshgrid([0, 8, 16, 24] * d, [0, 2, 4] * d);
            txPos = [txX(:).'; zeros(1, 12); txZ(:).'];
            
            arrayObj = chirp.sim.MimoArray(txPos, rxPos);
        end
    end
end
