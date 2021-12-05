classdef SwerlingModel
    % SWERLINGMODEL Fluctuation statistics for radar target cross sections.
    %
    % Generates stochastic radar cross section (RCS) realisations according to
    % the classical Swerling fluctuation cases:
    %   Swerling 0 / Constant: Non-fluctuating target
    %   Swerling 1: Rayleigh amplitude, slow scan-to-scan fluctuation
    %   Swerling 2: Rayleigh amplitude, fast pulse-to-pulse fluctuation
    %   Swerling 3: Chi-square with 4 degrees of freedom, slow fluctuation
    %   Swerling 4: Chi-square with 4 degrees of freedom, fast fluctuation
    
    methods (Static)
        function rcsSamples = sample(swerlingCase, meanRcs, numChirps)
            % Generate RCS realisations across slow-time chirps in a frame.
            %
            % Inputs:
            %   swerlingCase - Integer index 0, 1, 2, 3, or 4
            %   meanRcs      - Mean physical cross section in square meters
            %   numChirps    - Number of slow-time pulses/chirps in frame
            %
            % Output:
            %   rcsSamples   - 1 x numChirps positive double vector
            
            if meanRcs <= 0
                error('SwerlingModel:InvalidMeanRcs', 'Mean RCS must be strictly positive.');
            end
            
            switch swerlingCase
                case 0
                    % Non-fluctuating target (constant RCS)
                    rcsSamples = meanRcs * ones(1, numChirps);
                    
                case 1
                    % Swerling 1: Exponential power distribution, constant per frame
                    % p(sigma) = (1 / sigma_avg) * exp(-sigma / sigma_avg)
                    singleVal = -meanRcs * log(max(rand(), 1e-12));
                    rcsSamples = singleVal * ones(1, numChirps);
                    
                case 2
                    % Swerling 2: Independent exponential distribution per chirp
                    u = max(rand(1, numChirps), 1e-12);
                    rcsSamples = -meanRcs * log(u);
                    
                case 3
                    % Swerling 3: Chi-square (4 DOF), constant per frame
                    % Sum of two exponential variables each with mean meanRcs / 2
                    u1 = max(rand(), 1e-12);
                    u2 = max(rand(), 1e-12);
                    singleVal = -0.5 * meanRcs * (log(u1) + log(u2));
                    rcsSamples = singleVal * ones(1, numChirps);
                    
                case 4
                    % Swerling 4: Chi-square (4 DOF), independent per chirp
                    u1 = max(rand(1, numChirps), 1e-12);
                    u2 = max(rand(1, numChirps), 1e-12);
                    rcsSamples = -0.5 * meanRcs * (log(u1) + log(u2));
                    
                otherwise
                    error('SwerlingModel:UnknownCase', ...
                        'Swerling case %d not recognized. Must be 0, 1, 2, 3, or 4.', swerlingCase);
            end
        end
    end
end
