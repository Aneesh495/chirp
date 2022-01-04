classdef CaCfar < chirp.cfar.CfarDetector2D
    % CACFAR Specialized Cell-Averaging Constant False Alarm Rate detector.
    %
    % Assumes homogeneous background clutter and averages all reference training cells.
    
    methods
        function obj = CaCfar(pfa, trR, gdR, trV, gdV)
            % Construct a CaCfar detector.
            obj@chirp.cfar.CfarDetector2D();
            obj.Mode = 'CA';
            if nargin >= 1 && ~isempty(pfa), obj.Pfa = pfa; end
            if nargin >= 2 && ~isempty(trR), obj.TrainRange = trR; end
            if nargin >= 3 && ~isempty(gdR), obj.GuardRange = gdR; end
            if nargin >= 4 && ~isempty(trV), obj.TrainVelocity = trV; end
            if nargin >= 5 && ~isempty(gdV), obj.GuardVelocity = gdV; end
        end
    end
end
