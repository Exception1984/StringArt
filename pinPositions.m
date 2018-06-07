function [pinPos, angles] = pinPositions(numPins)
    
    if nargin == 0
        numPins = 512;
    end
    
    maxAngle = 2 * pi /numPins;
    minAngle = 0.95 * 2 * pi /numPins;
    range = maxAngle - minAngle;
    
    rng(0);
    
    % Accumulative Method
    angles = range * rand(numPins + 1, 1) + minAngle;
    angles = cumsum(angles);
    firstHookAngle = 2 * pi + angles(1);
    angles = angles .* (firstHookAngle / angles(end));
    angles = angles(1 : end - 1);
    
    pinPos = [cos(angles) sin(angles)];
end