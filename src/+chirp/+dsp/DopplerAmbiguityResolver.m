classdef DopplerAmbiguityResolver
    % DOPPLERAMBIGUITYRESOLVER Velocity de-aliasing via dual-PRF staggered chirps.
    %
    % Resolves radial velocity ambiguities when target velocity exceeds the Nyquist
    % unambiguous limit v_max = lambda / (4 * Tc) by combining two staggered PRF channels.
    
    properties
        Lambda                        % Operating carrier wavelength (m)
        Prf1                          % Pulse repetition frequency 1 (Hz)
        Prf2                          % Pulse repetition frequency 2 (Hz)
    end
    
    properties (Dependent)
        Vmax1                         % Max unambiguous velocity PRF 1 (m/s)
        Vmax2                         % Max unambiguous velocity PRF 2 (m/s)
        ExpandedVmax                  % De-aliased extended maximum velocity (m/s)
    end
    
    methods
        function obj = DopplerAmbiguityResolver(lambdaMeters, prf1Hz, prf2Hz)
            % Construct a dual-PRF Doppler ambiguity resolver.
            if nargin >= 3
                obj.Lambda = lambdaMeters;
                obj.Prf1 = prf1Hz;
                obj.Prf2 = prf2Hz;
            end
        end
        
        function v = get.Vmax1(obj)
            v = (obj.Lambda * obj.Prf1) / 4.0;
        end
        
        function v = get.Vmax2(obj)
            v = (obj.Lambda * obj.Prf2) / 4.0;
        end
        
        function vExp = get.ExpandedVmax(obj)
            % Expanded unambiguous velocity: V_exp = (Vmax1 * Vmax2) / abs(Vmax1 - Vmax2)
            dv = abs(obj.Vmax1 - obj.Vmax2);
            if dv < 1e-6
                vExp = obj.Vmax1;
            else
                vExp = (obj.Vmax1 * obj.Vmax2) / dv;
            end
        end
        
        function trueVelocity = resolveChineseRemainder(obj, measuredV1, measuredV2, maxFolds)
            % Resolve true velocity using Chinese Remainder Theorem integer search.
            %
            % Inputs:
            %   measuredV1 - Aliased velocity measured in PRF 1 (m/s)
            %   measuredV2 - Aliased velocity measured in PRF 2 (m/s)
            %   maxFolds   - Maximum integer ambiguity hypothesis search bounds (default: 5)
            %
            % Returns:
            %   trueVelocity - De-aliased true kinematic radial velocity (m/s)
            
            if nargin < 4 || isempty(maxFolds)
                maxFolds = 5;
            end
            
            v1_span = 2.0 * obj.Vmax1;
            v2_span = 2.0 * obj.Vmax2;
            
            bestDiff = Inf;
            trueVelocity = measuredV1;
            
            % Search integer ambiguity hypothesis pairs (k1, k2)
            for k1 = -maxFolds:maxFolds
                vCand1 = measuredV1 + k1 * v1_span;
                for k2 = -maxFolds:maxFolds
                    vCand2 = measuredV2 + k2 * v2_span;
                    
                    diffVal = abs(vCand1 - vCand2);
                    if diffVal < bestDiff
                        bestDiff = diffVal;
                        trueVelocity = 0.5 * (vCand1 + vCand2);
                    end
                end
            end
        end
    end
end
