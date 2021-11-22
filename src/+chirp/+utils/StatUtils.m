classdef StatUtils
    % STATUTILS Statistical distributions and detection probability functions for radar.
    %
    % Provides probability density and cumulative distribution functions for
    % Rayleigh (Swerling I/II), Chi-square (Swerling III/IV), Weibull, K-distribution
    % clutter, and Marcum Q detection integrals.
    
    methods (Static)
        function p = rayleigh_pdf(x, sigma)
            % Rayleigh probability density function: (x / sigma^2) * exp(-x^2 / (2*sigma^2))
            if nargin < 2, sigma = 1.0; end
            p = zeros(size(x));
            mask = x >= 0;
            xm = x(mask);
            p(mask) = (xm ./ (sigma^2)) .* exp(-(xm.^2) ./ (2.0 * sigma^2));
        end
        
        function c = rayleigh_cdf(x, sigma)
            % Rayleigh cumulative distribution function: 1 - exp(-x^2 / (2*sigma^2))
            if nargin < 2, sigma = 1.0; end
            c = zeros(size(x));
            mask = x >= 0;
            xm = x(mask);
            c(mask) = 1.0 - exp(-(xm.^2) ./ (2.0 * sigma^2));
        end
        
        function p = weibull_pdf(x, a, b)
            % Weibull PDF: (b/a) * (x/a)^(b-1) * exp(-(x/a)^b)
            %   a - Scale parameter (a > 0)
            %   b - Shape parameter (b > 0)
            p = zeros(size(x));
            mask = x > 0;
            xm = x(mask);
            p(mask) = (b / a) .* ((xm ./ a).^(b - 1.0)) .* exp(-((xm ./ a).^b));
        end
        
        function p = chisquare_pdf(x, k)
            % Chi-square probability density function with k degrees of freedom.
            %   p(x; k) = (1 / (2^(k/2) * Gamma(k/2))) * x^(k/2 - 1) * exp(-x/2)
            p = zeros(size(x));
            mask = x > 0;
            xm = x(mask);
            c = 1.0 / (2.0^(k / 2.0) * gamma(k / 2.0));
            p(mask) = c .* (xm.^(k / 2.0 - 1.0)) .* exp(-xm ./ 2.0);
        end
        
        function q = q_func(x)
            % Gaussian Q-function: Q(x) = 0.5 * erfc(x / sqrt(2))
            q = 0.5 * erfc(x ./ sqrt(2.0));
        end
        
        function prob = marcum_q(a, b, m)
            % Generalized Marcum Q-function Q_m(a, b) approximation.
            %
            % Computes probability that normalized envelope of sinusoid plus narrow-band
            % Gaussian noise exceeds threshold b, with noncentrality parameter a.
            % m defaults to 1 for first-order Marcum Q.
            
            if nargin < 3, m = 1; end
            
            % High-precision series summation approximation
            if a == 0
                % Central chi-square upper tail
                prob = exp(-0.5 * b.^2);
                return;
            end
            
            numTerms = 35;
            term = exp(-0.5 * (a.^2 + b.^2));
            sumSeries = 0.0;
            
            % Modified Bessel function series expansion
            for k = (1 - m):numTerms
                sumSeries = sumSeries + ((a ./ b).^k) .* besseli(k, a .* b);
            end
            
            prob = term .* sumSeries;
            prob = max(0.0, min(1.0, prob));
        end
    end
end
