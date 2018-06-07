function [edgeCodes, HookToEdge] = edgeCodes(numPins)
    hookIndices = (1 : 1 : numPins)';
    edgeCodes = nchoosek(hookIndices, 2);
    edgeCodes = [repelem(edgeCodes(:, 1), 4) repelem(edgeCodes(:, 2), 4)];
    
    % Consider four connections between each hook pair
    numConn = 4;

    numEdges = size(edgeCodes, 1);
    
    connectionType = repmat((1 : numConn)', numEdges, 1);
    
    hookToEdgeFinder = HookToEdgeFinder(edgeCodes, connectionType);
	HookToEdge = arrayfun(@(x) compute(hookToEdgeFinder, x), hookIndices, 'UniformOutput', false);
end