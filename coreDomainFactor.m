function fac = coreDomainFactor(minAngle)

    if nargin < 1
        minAngle = 7 * pi / 8;
    end

    fac = max(0.0, min(1.0, sin(0.5 * (pi - minAngle))));
end