classdef RadarSimulator
    % RADARSIMULATOR Master physical simulation engine for multi-channel FMCW radar cubes.
    %
    % Synthesizes raw de-chirped baseband I/Q radar cubes [N_adc x N_chirps x N_virtual]
    % incorporating MIMO array steering, Doppler shifts, Swerling fluctuations,
    % atmospheric attenuation, and thermal AWGN noise.
    
    properties
        Config                % chirp.sim.RadarConfig instance
        Chirp                 % chirp.waveform.FmcwChirp instance
        Array                 % chirp.sim.MimoArray instance
        EnableNoise = true;   % Boolean flag to add thermal AWGN noise
        EnableQuantization = false; % Boolean flag to apply ADC quantization
    end
    
    methods
        function obj = RadarSimulator(configObj, chirpObj, arrayObj)
            % Construct a RadarSimulator instance with configuration objects.
            if nargin < 1 || isempty(configObj)
                obj.Config = chirp.sim.RadarConfig();
            else
                obj.Config = configObj;
            end
            
            if nargin < 2 || isempty(chirpObj)
                obj.Chirp = chirp.waveform.WaveformGenerator.automotiveLRR();
            else
                obj.Chirp = chirpObj;
            end
            
            if nargin < 3 || isempty(arrayObj)
                obj.Array = chirp.sim.MimoArray.automotive3Tx4Rx(obj.Chirp.Wavelength);
            else
                obj.Array = arrayObj;
            end
        end
        
        function [radarCube, timeGrid] = simulateFrame(obj, targetList, clutterPatches)
            % Synthesize a complete frame of raw baseband radar cube data.
            %
            % Inputs:
            %   targetList     - Array of chirp.sim.TargetModel objects
            %   clutterPatches - (Optional) Array of clutter TargetModel objects
            %
            % Output:
            %   radarCube      - [N_adc x N_chirps x N_virtual] complex double array
            %   timeGrid       - Struct with fastTime and slowTime coordinate vectors
            
            nAdc = obj.Chirp.NumAdcSamples;
            nChirps = obj.Config.NumChirps;
            nVirt = obj.Array.NumVirtual;
            wavelength = obj.Chirp.Wavelength;
            tChirp = obj.Chirp.TotalChirpPeriod;
            
            radarCube = zeros(nAdc, nChirps, nVirt);
            
            allTargets = targetList;
            if nargin >= 3 && ~isempty(clutterPatches)
                allTargets = [allTargets, clutterPatches];
            end
            
            numTgts = length(allTargets);
            
            % Pre-sample Swerling RCS across all chirps for each target
            rcsMatrix = zeros(numTgts, nChirps);
            for iTgt = 1:numTgts
                rcsMatrix(iTgt, :) = chirp.sim.SwerlingModel.sample(...
                    allTargets(iTgt).SwerlingType, ...
                    allTargets(iTgt).Rcs, ...
                    nChirps ...
                );
            end
            
            fastTime = obj.Chirp.fastTimeVector();
            slowTime = (0:(nChirps - 1)) * tChirp;
            
            % Synthesize signal across chirps
            for iChirp = 1:nChirps
                tSlow = slowTime(iChirp);
                
                % Chirp-level 2D slice [N_adc x N_virtual]
                chirpSlice = zeros(nAdc, nVirt);
                
                for iTgt = 1:numTgts
                    tgt = allTargets(iTgt);
                    
                    % Propagate target kinematic position to current chirp epoch
                    tgtCurrent = tgt.propagate(tSlow);
                    r = tgtCurrent.Range;
                    vr = tgtCurrent.RadialVelocity;
                    az = tgtCurrent.AzimuthDeg;
                    el = tgtCurrent.ElevationDeg;
                    
                    if r < 0.5 || r > obj.Chirp.MaxUnambiguousRange
                        continue;
                    end
                    
                    % Received power and amplitude
                    rcsVal = rcsMatrix(iTgt, iChirp);
                    prx = obj.Config.receivedPower(r, rcsVal, wavelength);
                    
                    % Atmospheric two-way loss factor
                    lossAtm = chirp.sim.AtmosphericLoss.twoWayLoss(r, obj.Chirp.StartFrequency);
                    prx = prx * lossAtm;
                    
                    voltageAmp = sqrt(max(prx, 0.0));
                    
                    % 1D fast-time beat signal vector [1 x N_adc]
                    beatSig = obj.Chirp.dechirpTarget(r, vr, voltageAmp);
                    
                    % Virtual MIMO spatial steering vector [N_virt x 1]
                    sv = obj.Array.steeringVector(az, el, wavelength);
                    
                    % Accumulate outer product: [N_adc x 1] * [1 x N_virt]
                    chirpSlice = chirpSlice + (beatSig(:) * sv.');
                end
                
                radarCube(:, iChirp, :) = reshape(chirpSlice, [nAdc, 1, nVirt]);
            end
            
            % Add calibrated receiver thermal AWGN noise
            if obj.EnableNoise
                noisePower = obj.Config.NoiseFloorWatts;
                noiseSigma = sqrt(noisePower / 2.0);
                noise = noiseSigma * (randn(size(radarCube)) + 1j * randn(size(radarCube)));
                radarCube = radarCube + noise;
            end
            
            % ADC Quantization
            if obj.EnableQuantization
                maxVal = max(abs(radarCube(:)));
                if maxVal > 0
                    qLevels = 2^(obj.Config.AdcBits - 1) - 1;
                    scale = qLevels / maxVal;
                    radarCube = (round(real(radarCube) * scale) + 1j * round(imag(radarCube) * scale)) / scale;
                end
            end
            
            timeGrid = struct('fastTime', fastTime, 'slowTime', slowTime);
        end
    end
end
