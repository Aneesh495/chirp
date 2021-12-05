classdef RadarConfig
    % RADARCONFIG System configuration and link budget parameters for FMCW radar.
    %
    % Encapsulates transmit power, receiver noise figure, antenna gains,
    % loss factors, and ADC quantization parameters.
    
    properties
        TxPowerDbm = 12.0;             % Transmit power in dBm (e.g. 12 dBm ~ 16 mW)
        TxAntennaGainDbi = 10.0;       % Transmit antenna gain in dBi
        RxAntennaGainDbi = 10.0;       % Receive antenna gain in dBi
        NoiseFigureDb = 12.0;          % Receiver noise figure in dB
        SystemLossDb = 3.0;            % Implementation and RF interconnect loss in dB
        AdcBits = 16;                  % ADC resolution in bits (e.g. 12, 14, or 16)
        AdcFullScaleVoltage = 1.0;     % ADC full-scale voltage peak (V)
        TemperatureKelvin = 290.0;     % Operating physical temperature (K)
        NumChirps = 64;                % Slow-time chirps per frame
        FramePeriod = 50.0e-3;         % Frame repetition period (s, 20 Hz frame rate)
    end
    
    properties (Dependent)
        TxPowerWatts          % Transmit power in Watts
        NoiseFloorWatts       % Total receiver thermal noise floor power
        NoiseFloorDbm         % Thermal noise floor in dBm
    end
    
    methods
        function obj = RadarConfig(varargin)
            % Construct a RadarConfig instance with optional name-value overrides.
            if nargin > 0
                for k = 1:2:nargin
                    propName = varargin{k};
                    propVal = varargin{k+1};
                    if isprop(obj, propName)
                        obj.(propName) = propVal;
                    else
                        error('RadarConfig:UnknownProperty', 'Property %s not recognized.', propName);
                    end
                end
            end
        end
        
        function p = get.TxPowerWatts(obj)
            p = chirp.utils.PhysicalConstants.dbm2watts(obj.TxPowerDbm);
        end
        
        function nf = get.NoiseFloorWatts(obj)
            % P_noise = k * T * B * F_n
            % Default evaluation with reference 25 MHz IF bandwidth
            k = chirp.utils.PhysicalConstants.BOLTZMANN;
            fn = chirp.utils.PhysicalConstants.db2pow(obj.NoiseFigureDb);
            bandwidth = 25.0e6;
            nf = k * obj.TemperatureKelvin * bandwidth * fn;
        end
        
        function nfd = get.NoiseFloorDbm(obj)
            nfd = chirp.utils.PhysicalConstants.watts2dbm(obj.NoiseFloorWatts);
        end
        
        function prxWatts = receivedPower(obj, rangeMeters, rcsM2, wavelengthMeters)
            % Compute received power using the classic radar range equation:
            %   Pr = (Pt * Gt * Gr * lambda^2 * sigma) / ((4*pi)^3 * R^4 * L)
            
            c = chirp.utils.PhysicalConstants.C;
            pt = obj.TxPowerWatts;
            gt = chirp.utils.PhysicalConstants.db2pow(obj.TxAntennaGainDbi);
            gr = chirp.utils.PhysicalConstants.db2pow(obj.RxAntennaGainDbi);
            loss = chirp.utils.PhysicalConstants.db2pow(obj.SystemLossDb);
            
            numerator = pt * gt * gr * (wavelengthMeters^2) * rcsM2;
            denominator = ((4.0 * pi)^3) * (rangeMeters.^4) * loss;
            
            prxWatts = numerator ./ max(denominator, eps);
        end
    end
end
