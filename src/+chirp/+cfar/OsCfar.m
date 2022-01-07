classdef OsCfar < chirp.cfar.CfarDetector2D
    % OSCFAR Specialized Ordered-Statistic Constant False Alarm Rate detector.
    %
    % Sorts reference training cells and selects the k-th percentile value
    % to prevent target masking in dense multi-target and inhomogeneous clutter scenes.
    
    methods
        function obj = OsCfar(pfa, trR, gdR, trV, gdV, rankK)
            % Construct an OsCfar detector.
            obj@chirp.cfar.CfarDetector2D();
            obj.Mode = 'OS';
            if nargin >= 1 && ~isempty(pfa), obj.Pfa = pfa; end
            if nargin >= 2 && ~isempty(trR), obj.TrainRange = trR; end
            if nargin >= 3 && ~isempty(gdR), obj.GuardRange = gdR; end
            if nargin >= 4 && ~isempty(trV), obj.TrainVelocity = trV; end
            if nargin >= 5 && ~isempty(gdV), obj.GuardVelocity = gdV; end
            if nargin >= 6 && ~isempty(rankK), obj.RankK = rankK; end
        end
    end
end
