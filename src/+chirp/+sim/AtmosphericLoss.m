classdef AtmosphericLoss
    % ATMOSPHERICLOSS Electromagnetic attenuation through atmospheric gases.
    %
    % Computes specific attenuation in dB/km and two-way path loss for mmWave radar
    % based on simplified ITU-R P.676 resonance lines for oxygen and water vapor.
    
    methods (Static)
        function gammaDbPerKm = specificAttenuation(freqHz, pressureHpa, tempC, humidityGPerM3)
            % Compute specific attenuation in dB/km at standard sea level.
            %
            % Inputs:
            %   freqHz         - Carrier frequency in Hz (e.g. 77e9)
            %   pressureHpa    - Dry atmospheric pressure in hPa (default: 1013.25)
            %   tempC          - Temperature in Celsius (default: 15.0 C)
            %   humidityGPerM3 - Water vapor concentration in g/m^3 (default: 7.5 g/m^3)
            
            if nargin < 2 || isempty(pressureHpa), pressureHpa = 1013.25; end
            if nargin < 3 || isempty(tempC), tempC = 15.0; end
            if nargin < 4 || isempty(humidityGPerM3), humidityGPerM3 = 7.5; end
            
            fGhz = freqHz / 1.0e9;
            theta = 300.0 / (273.15 + tempC);
            
            % Oxygen attenuation: strong resonance line near 60 GHz
            fO2 = 60.0;
            widthO2 = 4.0;
            gammaO2 = 0.005 * (fGhz^2) * (theta^2) * (pressureHpa / 1013.25) + ...
                      15.0 * exp(-((fGhz - fO2) / widthO2)^2);
            
            % Water vapor attenuation: resonance lines near 22.2 GHz and 183 GHz
            fH2O = 22.235;
            widthH2O = 3.0;
            gammaH2O = 0.0002 * (fGhz^2) * (humidityGPerM3) + ...
                       0.2 * (humidityGPerM3 / 7.5) * exp(-((fGhz - fH2O) / widthH2O)^2);
            
            % Total specific attenuation (one-way, dB/km)
            gammaDbPerKm = gammaO2 + gammaH2O;
        end
        
        function lossFactor = twoWayLoss(rangeMeters, freqHz)
            % Compute linear two-way round-trip atmospheric transmission power loss factor.
            % Returns a fraction <= 1.0.
            
            gamma = chirp.sim.AtmosphericLoss.specificAttenuation(freqHz);
            rangeKm = rangeMeters / 1000.0;
            twoWayLossDb = 2.0 * gamma * rangeKm;
            
            lossFactor = 10.0^(-twoWayLossDb / 10.0);
        end
    end
end
