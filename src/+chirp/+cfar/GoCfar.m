classdef GoCfar < chirp.cfar.CfarDetector2D
    % GOCFAR Greatest-Of Constant False Alarm Rate detector.
    %
    % Takes the maximum of leading and lagging training half-windows to maintain
    % low false alarm rates across sharp clutter transition boundaries.
    
    methods
        function obj = GoCfar(pfa, trR, gdR, trV, gdV)
            % Construct a GoCfar detector.
            obj@chirp.cfar.CfarDetector2D();
            obj.Mode = 'GO';
            if nargin >= 1 && ~isempty(pfa), obj.Pfa = pfa; end
            if nargin >= 2 && ~isempty(trR), obj.TrainRange = trR; end
            if nargin >= 3 && ~isempty(gdR), obj.GuardRange = gdR; end
            if nargin >= 4 && ~isempty(trV), obj.TrainVelocity = trV; end
            if nargin >= 5 && ~isempty(gdV), obj.GuardVelocity = gdV; end
        end
    end
end
