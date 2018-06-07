function [e, nout] = computeEdgeInPositiveDomain(x, y, n, domainWidth)

e = computePixelCodeInPositiveDomain(x, y, domainWidth);
nout = n .* ones(size(x,1), 1);

end