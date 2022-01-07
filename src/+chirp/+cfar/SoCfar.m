classdef SoCfar < chirp.cfar.CfarDetector2D
    % SOCFAR Smallest-Of Constant False Alarm Rate detector.
    %
    % Takes the minimum of leading and lagging training half-windows to resolve
    % closely spaced adjacent targets without mutual target masking.
    
    methods
        function obj = SoCfar(pfa, trR, gdR, trV, gdV)
            % Construct a SoCfar detector.
            obj@chirp.cfar.CfarDetector2D();
            obj.Mode = 'SO';
            if nargin >= 1 && ~isempty(pfa), obj.Pfa = pfa; end
            if nargin >= 2 && ~isempty(trR), obj.TrainRange = trR; end
            if nargin >= 3 && ~isempty(gdR), obj.GuardRange = gdR; end
            if nargin >= 4 && ~isempty(trV), obj.TrainVelocity = trV; end
            if nargin >= 5 && ~isempty(gdV), obj.GuardVelocity = gdV; end
        end
    end
end
