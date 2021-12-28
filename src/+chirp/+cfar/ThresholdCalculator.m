classdef ThresholdCalculator
    % THRESHOLDCALCULATOR Analytical threshold scaling factors for CFAR detection.
    %
    % Computes alpha multipliers that guarantee specified false alarm probability P_fa
    % for Cell-Averaging (CA), Ordered-Statistic (OS), Greatest-Of (GO), and Smallest-Of (SO) CFAR.
    
    methods (Static)
        function alpha = caCfarAlpha(pfa, numTrainCells)
            % Compute scaling factor alpha for Cell-Averaging CFAR.
            % Formula: alpha = N * (Pfa^(-1/N) - 1)
            
            if pfa <= 0 || pfa >= 1
                error('ThresholdCalculator:InvalidPfa', 'Pfa must be in (0, 1).');
            end
            if numTrainCells <= 0
                error('ThresholdCalculator:InvalidTrainCells', 'numTrainCells must be positive.');
            end
            
            N = numTrainCells;
            alpha = N * (pfa^(-1.0 / N) - 1.0);
        end
        
        function alpha = osCfarAlpha(pfa, numTrainCells, kRank)
            % Compute scaling factor alpha for Ordered-Statistic CFAR via bisection root finding.
            %
            % Pfa formula:
            %   Pfa = k * nchoosek(N, k) * Gamma(k) * Gamma(N - k + 1 + alpha) / Gamma(N + 1 + alpha)
            %       = k! * nchoosek(N, k) / prod_{j=0}^{k-1} (N - j + alpha)
            
            if nargin < 3 || isempty(kRank)
                kRank = round(0.75 * numTrainCells);
            end
            
            N = numTrainCells;
            k = max(1, min(N, round(kRank)));
            
            % Objective function: log(Pfa_calc(alpha)) - log(pfa) = 0
            % log(Pfa) = sum_{j=1}^k log(N - j + 1) - sum_{j=0}^{k-1} log(N - j + alpha)
            logPfaTarget = log(pfa);
            logNum = sum(log((N - k + 1):N));
            
            f = @(a) logNum - sum(log((N - (0:(k - 1))) + a)) - logPfaTarget;
            
            % Bracket root [0.01, 1000]
            lowA = 0.01;
            highA = 200.0;
            while f(highA) > 0
                highA = highA * 2.0;
            end
            
            % Bisection search
            for iter = 1:50
                midA = 0.5 * (lowA + highA);
                val = f(midA);
                if abs(val) < 1e-9 || (highA - lowA) < 1e-6
                    alpha = midA;
                    return;
                end
                if val > 0
                    lowA = midA;
                else
                    highA = midA;
                end
            end
            alpha = midA;
        end
        
        function alpha = goCfarAlpha(pfa, numTrainCells)
            % Approximate scaling factor alpha for Greatest-Of CFAR.
            % Leading and lagging halves have M = N / 2 training cells each.
            
            M = max(1, round(numTrainCells / 2.0));
            % Pfa = 2 * (1 + alpha/M)^(-M) - (1 + 2*alpha/M)^(-M)
            f = @(a) 2.0 * ((1.0 + a / M)^(-M)) - ((1.0 + 2.0 * a / M)^(-M)) - pfa;
            
            lowA = 0.1;
            highA = 200.0;
            while f(highA) > 0
                highA = highA * 2.0;
            end
            
            for iter = 1:50
                midA = 0.5 * (lowA + highA);
                val = f(midA);
                if abs(val) < 1e-9 || (highA - lowA) < 1e-6
                    alpha = midA;
                    return;
                end
                if val > 0
                    lowA = midA;
                else
                    highA = midA;
                end
            end
            alpha = midA;
        end
        
        function alpha = soCfarAlpha(pfa, numTrainCells)
            % Approximate scaling factor alpha for Smallest-Of CFAR.
            M = max(1, round(numTrainCells / 2.0));
            % Pfa = 2 * (1 + alpha/M)^(-M) - (1 + alpha/M)^(-2*M)
            f = @(a) 2.0 * ((1.0 + a / M)^(-M)) - ((1.0 + a / M)^(-2.0 * M)) - pfa;
            
            lowA = 0.1;
            highA = 200.0;
            while f(highA) > 0
                highA = highA * 2.0;
            end
            
            for iter = 1:50
                midA = 0.5 * (lowA + highA);
                val = f(midA);
                if abs(val) < 1e-9 || (highA - lowA) < 1e-6
                    alpha = midA;
                    return;
                end
                if val > 0
                    lowA = midA;
                else
                    highA = midA;
                end
            end
            alpha = midA;
        end
    end
end
