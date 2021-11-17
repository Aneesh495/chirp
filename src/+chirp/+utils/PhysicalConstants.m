classdef PhysicalConstants
    % PHYSICALCONSTANTS Fundamental physical constants and electromagnetic parameters.
    %
    % This class provides SI constants, reference temperatures, and standard
    % conversion utilities for radar systems and electromagnetic wave propagation.
    %
    % Example:
    %   c = chirp.utils.PhysicalConstants.C;
    %   lambda = chirp.utils.PhysicalConstants.freq2wavelength(77e9);
    
    properties (Constant)
        % Speed of light in vacuum (m/s)
        C = 299792458.0;
        
        % Boltzmann constant (J/K)
        BOLTZMANN = 1.380649e-23;
        
        % Standard reference noise temperature (Kelvin)
        T0 = 290.0;
        
        % Permittivity of free space (F/m)
        EPSILON_0 = 8.8541878128e-12;
        
        % Permeability of free space (H/m)
        MU_0 = 1.25663706212e-6;
        
        % Characteristic impedance of free space (Ohms)
        ETA_0 = 376.730313668;
    end
    
    methods (Static)
        function lambda = freq2wavelength(freqHz)
            % Convert carrier frequency in Hertz to wavelength in meters.
            if any(freqHz <= 0)
                error('PhysicalConstants:InvalidFrequency', 'Frequency must be strictly positive.');
            end
            lambda = chirp.utils.PhysicalConstants.C ./ freqHz;
        end
        
        function freqHz = wavelength2freq(lambdaMeters)
            % Convert wavelength in meters to carrier frequency in Hertz.
            if any(lambdaMeters <= 0)
                error('PhysicalConstants:InvalidWavelength', 'Wavelength must be strictly positive.');
            end
            freqHz = chirp.utils.PhysicalConstants.C ./ lambdaMeters;
        end
        
        function lin = db2pow(db)
            % Convert decibels to power linear scale: 10^(db/10)
            lin = 10.^(db ./ 10.0);
        end
        
        function db = pow2db(lin)
            % Convert power linear scale to decibels: 10*log10(lin)
            linSafe = max(lin, eps);
            db = 10.0 .* log10(linSafe);
        end
        
        function lin = db2mag(db)
            % Convert decibels to magnitude linear scale: 10^(db/20)
            lin = 10.^(db ./ 20.0);
        end
        
        function db = mag2db(lin)
            % Convert magnitude linear scale to decibels: 20*log10(lin)
            linSafe = max(abs(lin), eps);
            db = 20.0 .* log10(linSafe);
        end
        
        function watts = dbm2watts(dbm)
            % Convert power in dBm to Watts: 10^((dbm - 30)/10)
            watts = 10.^((dbm - 30.0) ./ 10.0);
        end
        
        function dbm = watts2dbm(watts)
            % Convert power in Watts to dBm: 10*log10(watts) + 30
            wattsSafe = max(watts, eps);
            dbm = 10.0 .* log10(wattsSafe) + 30.0;
        end
        
        function [r, az, el] = cart2sphDeg(x, y, z)
            % Convert Cartesian coordinates (x, y, z) to range, azimuth (deg), elevation (deg).
            % Azimuth is measured from +y axis towards +x (sensor boresight along +y),
            % or standard radar convention where x is cross-range and y is range.
            r = sqrt(x.^2 + y.^2 + z.^2);
            az = atan2(x, y) .* (180.0 / pi);
            el = atan2(z, sqrt(x.^2 + y.^2)) .* (180.0 / pi);
        end
        
        function [x, y, z] = sph2cartDeg(r, az, el)
            % Convert range, azimuth (deg), elevation (deg) to Cartesian (x, y, z).
            azRad = az .* (pi / 180.0);
            elRad = el .* (pi / 180.0);
            x = r .* cos(elRad) .* sin(azRad);
            y = r .* cos(elRad) .* cos(azRad);
            z = r .* sin(elRad);
        end
    end
end
