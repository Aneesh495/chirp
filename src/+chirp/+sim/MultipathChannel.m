classdef MultipathChannel
    % MULTIPATHCHANNEL Ground bounce specular multipath reflection model.
    %
    % Simulates two-ray ground reflection producing constructive and destructive
    % interference fringes, as well as mirror ghost target returns.
    
    properties
        SurfaceType = 'asphalt';     % Surface: 'asphalt', 'concrete', or 'dry-ground'
        RoughnessMeters = 0.005;     % Surface RMS height roughness (m)
        SensorHeight = 0.8;          % Radar mount height above ground (m)
    end
    
    methods
        function obj = MultipathChannel(varargin)
            % Construct a MultipathChannel instance with optional overrides.
            if nargin > 0
                for k = 1:2:nargin
                    propName = varargin{k};
                    propVal = varargin{k+1};
                    if isprop(obj, propName)
                        obj.(propName) = propVal;
                    else
                        error('MultipathChannel:UnknownProperty', 'Property %s not found.', propName);
                    end
                end
            end
        end
        
        function [gammaH, gammaV] = reflectionCoefficients(obj, grazingAngleRad, freqHz)
            % Compute horizontal (TE) and vertical (TM) Fresnel reflection coefficients.
            
            switch lower(obj.SurfaceType)
                case 'asphalt'
                    epsR = 4.5;
                    sigmaCond = 0.02;
                case 'concrete'
                    epsR = 6.5;
                    sigmaCond = 0.05;
                case 'dry-ground'
                    epsR = 3.0;
                    sigmaCond = 0.01;
                otherwise
                    epsR = 4.0;
                    sigmaCond = 0.02;
            end
            
            % Complex relative permittivity: eps_c = eps_r - 1j * sigma / (2*pi*f*eps0)
            eps0 = chirp.utils.PhysicalConstants.EPSILON_0;
            epsC = epsR - 1j * sigmaCond / (2.0 * pi * freqHz * eps0);
            
            sinPsi = sin(grazingAngleRad);
            sqrtTerm = sqrt(epsC - cos(grazingAngleRad).^2);
            
            % Horizontal polarization
            gammaH = (sinPsi - sqrtTerm) ./ (sinPsi + sqrtTerm);
            
            % Vertical polarization
            gammaV = (epsC .* sinPsi - sqrtTerm) ./ (epsC .* sinPsi + sqrtTerm);
            
            % Ament surface roughness attenuation factor: rho = exp(-2 * (2*pi*sigma_h * sin(psi) / lambda)^2)
            c = chirp.utils.PhysicalConstants.C;
            lambda = c / freqHz;
            phaseRough = (4.0 * pi * obj.RoughnessMeters * sinPsi) / lambda;
            rho = exp(-0.5 * (phaseRough.^2));
            
            gammaH = gammaH * rho;
            gammaV = gammaV * rho;
        end
        
        function ghost = createGroundGhost(obj, directTarget)
            % Synthesize a ground-bounce image ghost target.
            % Mirror image reflection appears below ground with inverted z coordinate:
            %   z_ghost = -(2 * h_sensor + z_target)
            
            ghost = directTarget;
            h = obj.SensorHeight;
            
            % Target height above ground
            targetHeightAboveGround = h + directTarget.Position(3);
            
            % Mirror position relative to radar
            ghostPos = directTarget.Position;
            ghostPos(3) = -(h + targetHeightAboveGround);
            
            ghost.Position = ghostPos;
            % Attenuate RCS due to ground reflection loss
            ghost.Rcs = directTarget.Rcs * 0.25;
            ghost.TargetId = directTarget.TargetId + 50000;
        end
    end
end
