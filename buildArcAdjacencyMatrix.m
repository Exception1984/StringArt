function [hookPos, fabricatableMask, i, j, v, m, n] = buildArcAdjacencyMatrix(numPins, threadThickness, frameDiameter, pinSideLength)
    
    debugPlot = false;

    radiusMillimeters = 0.5 * frameDiameter;
    domainWidthPixels = round(frameDiameter / threadThickness);
    
    stringWidth = 2.0 / domainWidthPixels * radiusMillimeters;

    fprintf('Building Adjacency Matrix...\n');
       
    radiusPixels = 0.5 * domainWidthPixels;
    
	hooksUnknownDimension = pinPositions(numPins);
    
    r = sqrt(hooksUnknownDimension(:,1).* hooksUnknownDimension(:,1) + hooksUnknownDimension(:,2) .* hooksUnknownDimension(:,2));
	hooksPixels = radiusPixels * (hooksUnknownDimension ./ [r r]);
    hooksAngles = atan2(hooksPixels(:,2), hooksPixels(:,1));
    
    realToImageFactor = radiusPixels / radiusMillimeters;
    
    numPins = size(hooksUnknownDimension, 1);
    hookIndices = (1 : 1 : numPins)';
    edgeCodes = nchoosek(hookIndices, 2);
    
    hookArray = ObjectArray([repelem(realToImageFactor * pinSideLength, numPins, 1) hooksPixels hooksAngles]);
    
    if debugPlot
        % Plot Figure
        figure(1);
        hold on;
        daspect([1 1 1]);
    end
    
    numFabricatableEdges = 0;
    numEdgesAltogether = 0;
    
    % Consider four connections between each hook pair
    numConn = 4;
    
    % Compute Edges
	x1 = zeros(numConn * size(edgeCodes, 1), 1) - domainWidthPixels;
    y1 = x1;
    x2 = x1;
    y2 = x1;
    
    fabricatableMask = false(size(x1, 1), 1);
    
    if debugPlot
        % Plot Hooks
        for i = 1 : numPins
            hook = hookArray(i, 1).Value;
            hookPos = [hook.p1 hook.p2 hook.p3 hook.p4 hook.p1];
            hookPoints = [hook.p1 hook.p2 hook.p3 hook.p4];

            xHook = hookPos(1, :);
            yHook = hookPos(2, :);

            plot(xHook, yHook, 'blue');
            for k = 1 : 2
                w = circleDiameters(k);
                off = offs(k);
                for j = 1 : 4
                    rectangle('Position',[hookPoints(:, j)' - [off off] w w],'Curvature', [1 1], 'EdgeColor', 'blue');
                end
            end
        end
    end

    edgeIndex = 1;
    for i = 1 : size(edgeCodes, 1)
        aIndex = edgeCodes(i, 1);
        bIndex = edgeCodes(i, 2);
        
        % Hook A
        hookA = hookArray(aIndex, 1).Value;
        %hookAPos = [hookA.p1 hookA.p2 hookA.p3 hookA.p4 hookA.p1];
        %xA = hookAPos(1, :);
        %yA = hookAPos(2, :);
        %hookAPoints = [hookA.p1 hookA.p2 hookA.p3 hookA.p4];
        
        % Hook B
        hookB = hookArray(bIndex, 1).Value;
        %hookBPos = [hookB.p1 hookB.p2 hookB.p3 hookB.p4 hookB.p1];
        %hookBPoints = [hookB.p1 hookB.p2 hookB.p3 hookB.p4];
        
        % Compute Strings
        [stringAB, stringBA, diagAB, diagBA] = hookA.computeStrings(hookB, realToImageFactor * stringWidth);
        
        stringData = [stringAB' stringBA' diagAB' diagBA'];
        
        % Neighbors
        neighborAIndex = aIndex + 1;
        neighborBIndex = bIndex - 1;
        neighborCIndex = aIndex - 1;
        neighborDIndex = bIndex + 1;
        
        if (neighborAIndex > numPins) neighborAIndex = 1; end
        if (neighborBIndex == 0) neighborBIndex = numPins; end
        if (neighborCIndex == 0) neighborCIndex = numPins; end
        if (neighborDIndex > numPins) neighborDIndex = 1; end
        
        neighborA = hookArray(neighborAIndex, 1).Value;
        neighborB = hookArray(neighborBIndex, 1).Value;
        neighborC = hookArray(neighborCIndex, 1).Value;
        neighborD = hookArray(neighborDIndex, 1).Value;
        
        for p = 1 : 4
            p1 = stringData(:, 2 * (p - 1) + 1);
            p2 = stringData(:, 2 * (p - 1) + 2);

            % Test for fabricability

            if (~neighborA.intersectsString(p1, p2)) && (~neighborB.intersectsString(p1, p2)) && ...
              (~neighborC.intersectsString(p1, p2)) && (~neighborD.intersectsString(p1, p2))
                numFabricatableEdges = numFabricatableEdges + 1;
                
                x1(edgeIndex) = p1(1);
                y1(edgeIndex) = p1(2);
                x2(edgeIndex) = p2(1);
                y2(edgeIndex) = p2(2);
                
                fabricatableMask(edgeIndex) = true;
            end
            
            edgeIndex = edgeIndex + 1;
        end
        
        numEdgesAltogether = numEdgesAltogether + 4;
        
        if debugPlot
            if i == 53
                % Draw String Connections for the edge
                colors = {'cyan' 'magenta' 'green' 'yellow'};

                stringData = [stringAB' stringBA' diagAB' diagBA'];

                for p = 1 : 4
                    p1 = stringData(:, 2 * (p - 1) + 1);
                    p2 = stringData(:, 2 * (p - 1) + 2);

                    plot([p1(1) p2(1)]', [p1(2) p2(2)]', colors{p});

                end            
            end
        end
    end
        
    edgeCodes = [repelem(edgeCodes(:, 1), 4) repelem(edgeCodes(:, 2), 4)];

    hookPos = radiusMillimeters * (hooksUnknownDimension ./ [r r]);

    % Round x and y values
    x1 = round(x1 + 0.5 * (domainWidthPixels + 1));
    y1 = round(y1 + 0.5 * (domainWidthPixels + 1));
    x2 = round(x2 + 0.5 * (domainWidthPixels + 1));
    y2 = round(y2 + 0.5 * (domainWidthPixels + 1));

    numEdges = size(x1, 1);

    xOffset = zeros(numEdges, 1);
    yOffset = zeros(numEdges, 1);

    domainMin = ones(size(x1, 1), 1);
    domainMax = domainWidthPixels * ones(numEdges, 1);

    numLines = (1 : numEdges)';
    numLinesTotal = numEdges * ones(numEdges, 1);

    [x, y, v] = arrayfun(@drawLine, x1, y1, x2, y2, xOffset, yOffset, numLines, numLinesTotal, domainMin, domainMax, 'UniformOutput', false);

    num = num2cell((1 : 1 : size(x,1))');

    domainWidth = num2cell(domainWidthPixels * ones(numEdges, 1));
    [e, n] = cellfun(@computeEdgeInPositiveDomain, x, y, num, domainWidth, 'UniformOutput', false);

    i = double(cell2mat(e));
    j = cell2mat(n);
    v = cell2mat(v);
    m = domainWidthPixels * domainWidthPixels;
    n = size(edgeCodes, 1);

    fprintf('Done...\n');
end