function [x, eCodes, stringList, hookSequence, hookApproachLeftSide, outsideMask] = findConsecutivePath(x, matrixPath)
    tokens = strsplit(matrixPath, '_');
    numPins = str2double(tokens{end}(1:end-4));
    
    eCodes = edgeCodes(numPins);

    if ~isConnected(eCodes, x)
        error('Graph has to consist of just one connected component!');
    end
    
    % Column Indices for graphInfo matrix
    V1 = 2;
    V2 = 3;
    VALENCE = 4;
    
    numEdges = size(eCodes, 1);
    numHooks = max(max(eCodes));
    numHookConn = 4;    
    numPointToPointConnections = numEdges / numHookConn;
    
    nodeIds = (1 : 1 : numHooks)';
    
    graphInfo = [(1 : 1 : numEdges)' eCodes x];
    graphInfoNonEmpty = graphInfo(x > 0, :);

    aToBType = repmat([1 2 3 4]', numPointToPointConnections, 1);
    bToAType = repmat([2 1 3 4]', numPointToPointConnections, 1);
    
    usedEdgeCodes = eCodes(x == 1, :);
    hookConnWithType = [usedEdgeCodes(:) [aToBType(x == 1); bToAType(x == 1)]];
    
    hookConnTypeOne   = accumarray(hookConnWithType(hookConnWithType(:, 2) == 1, 1), 1, [numHooks, 1]);
    hookConnTypeTwo   = accumarray(hookConnWithType(hookConnWithType(:, 2) == 2, 1), 1, [numHooks, 1]);
    hookConnTypeThree = accumarray(hookConnWithType(hookConnWithType(:, 2) == 3, 1), 1, [numHooks, 1]);
    hookConnTypeFour  = accumarray(hookConnWithType(hookConnWithType(:, 2) == 4, 1), 1, [numHooks, 1]);
    
    hookConnFromLeft = hookConnTypeTwo + hookConnTypeFour;
    hookConnFromRight = hookConnTypeOne + hookConnTypeThree;
    
    sumValencePerNode = hookConnFromLeft + hookConnFromRight;
    hookConnDiff = hookConnFromLeft - hookConnFromRight;
    
    usedNodeIds = unique([graphInfoNonEmpty(:,V1); graphInfoNonEmpty(:,V2)]);
    numUsedNodeIds = size(usedNodeIds, 1);

    gTopRight = full(sparse(graphInfo(:, V1), graphInfo(:, V2), graphInfo(:, VALENCE), numHooks, numHooks));
    bBottomLeft = full(sparse(graphInfo(:, V2), graphInfo(:, V1), graphInfo(:, VALENCE), numHooks, numHooks));
    G = gTopRight + bBottomLeft;
    
    allConnWithIds = [nodeIds sumValencePerNode];
    oddConnections = allConnWithIds(hookConnDiff ~= 0, :);

    %rop = nchoosek(oddConnections(:,1), 2);
    %rop = rop(randperm(size(rop, 1)), :);
    
    oddPairs = nchoosek(oddConnections(:,1), 2);

    B = arrayfun(@(x,y) G(x, y), oddPairs(:, 1), oddPairs(:, 2));
    
    % Compute Edge IDs of odd pairs
    % rop == remainingOddPairs
    rop = [oddPairs B];
    
    % TODO: Find a better sorting of the edges than number of existing
    % strings
    rop = flipud(sortrows(rop,3));
    
    %c = [0; cumsum(numHookConn * (255 : -1 : 2)')];
    c = [0; cumsum(numHookConn * (numPins - 1 : -1 : 2)')];

    ropEdgeIdsMin = c(rop(:, 1)) + 1 + (rop(:, 2) - rop(:, 1) - 1) * numHookConn;
    ropEdgeIds = repmat(ropEdgeIdsMin, 1, 4) + repmat((0 : 1 : 3), size(rop, 1), 1);
    
	edgeFilter = true(size(ropEdgeIds));
    
    edgeRowFilter = any(edgeFilter, 2);
    
    rop = rop(edgeRowFilter, :);
    edgeFilter = edgeFilter(edgeRowFilter, :);
    ropEdgeIds = ropEdgeIds(edgeRowFilter, :);
    
    auxNeigh = zeros(numHooks, numHooks);
    auxNeigh = cat(3, auxNeigh, auxNeigh, auxNeigh, auxNeigh);
    
	added = 0;
    edgeId = 0;
    while size(rop, 1) > 1  % sum(abs(hookConnDiff)) > 2
        v1 = rop(1, 1);
        v2 = rop(1, 2);
        ind = [v1 v2]';
        edgeRange = ropEdgeIds(1, :);
        filter = edgeFilter(1, :);

        couldAdd = false;
        maxNum = min(abs(hookConnDiff(v1)), abs(hookConnDiff(v2)));

        if hookConnDiff(v1) < 0
            if hookConnDiff(v2) < 0
                % hookConnDiff(v1) < 0, hookConnDiff(v2) < 0
                % v1 wants 2 or 4, v2 wants 1 or 4 => 4
                if filter(4)
                    couldAdd = true;
                    edgeId = edgeRange(4);
                    
                    hookConnTypeFour(ind) = hookConnTypeFour(ind) + maxNum;
                    hookConnFromLeft(ind) = hookConnFromLeft(ind) + maxNum;
                    hookConnDiff(ind) = hookConnDiff(ind) + maxNum;
                    
                    auxNeigh(ind(1), ind(2), 4) = auxNeigh(ind(1), ind(2), 4) + maxNum;
                    auxNeigh(ind(2), ind(1), 4) = auxNeigh(ind(2), ind(1), 4) + maxNum;
                end
            else
                % hookConnDiff(v1) < 0, hookConnDiff(v2) > 0
                % v1 wants 2 or 4, v2 wants 2 or 3 => 2
                % Note: Conn. Type 2 appears as Type 1 to v2
                if filter(2)
                    couldAdd = true;
                    edgeId = edgeRange(2);
                    
                    hookConnTypeTwo(v1) = hookConnTypeTwo(v1) + maxNum;
                    hookConnTypeOne(v2) = hookConnTypeOne(v2) + maxNum;

                    hookConnFromLeft(v1) = hookConnFromLeft(v1) + maxNum;
                    hookConnFromRight(v2) = hookConnFromRight(v2) + maxNum;

                    hookConnDiff(v1) = hookConnDiff(v1) + maxNum;
                    hookConnDiff(v2) = hookConnDiff(v2) - maxNum;
                    
                    auxNeigh(ind(1), ind(2), 1) = auxNeigh(ind(1), ind(2), 1) + maxNum;
                    auxNeigh(ind(2), ind(1), 2) = auxNeigh(ind(2), ind(1), 2) + maxNum;
                end
            end
        else
            if hookConnDiff(v2) < 0
                % hookConnDiff(v1) > 0, hookConnDiff(v2) < 0
                % v1 wants 1 or 3, v2 wants 1 or 4 => 1
                % Note: Conn. Type 1 appears as Type 2 to v2
                if filter(1)
                    couldAdd = true;
                    edgeId = edgeRange(1);

                    hookConnTypeOne(v1) = hookConnTypeOne(v1) + maxNum;
                    hookConnTypeTwo(v2) = hookConnTypeTwo(v2) + maxNum;

                    hookConnFromRight(v1) = hookConnFromRight(v1) + maxNum;
                    hookConnFromLeft(v2) = hookConnFromLeft(v2) + maxNum;

                    hookConnDiff(v1) = hookConnDiff(v1) - maxNum;
                    hookConnDiff(v2) = hookConnDiff(v2) + maxNum;
                    
                    auxNeigh(ind(1), ind(2), 2) = auxNeigh(ind(1), ind(2), 2) + maxNum;
                    auxNeigh(ind(2), ind(1), 1) = auxNeigh(ind(2), ind(1), 1) + maxNum;
                end
            else
                % hookConnDiff(v1) > 0, hookConnDiff(v2) > 0
                % v1 wants 1 or 3, v2 wants 2 or 3 => 3
                if filter(3)
                    couldAdd = true;
                    edgeId = edgeRange(3);
                    
                    hookConnTypeThree(ind) = hookConnTypeThree(ind) + maxNum;
                    hookConnFromRight(ind) = hookConnFromRight(ind) + maxNum;
                    hookConnDiff(ind) = hookConnDiff(ind) - maxNum;
                    
                    auxNeigh(ind(1), ind(2), 3) = auxNeigh(ind(1), ind(2), 3) + maxNum;
                    auxNeigh(ind(2), ind(1), 3) = auxNeigh(ind(2), ind(1), 3) + maxNum;
                end
            end
        end
        
        if couldAdd
            newValence = G(v1, v2) + maxNum;
            G(v1, v2) = newValence;
            G(v2, v1) = newValence;
            
            x(edgeId) = x(edgeId) + maxNum;
            
            added = added + maxNum;
            
            rem = false(size(rop(:, 1)));
            
            if hookConnDiff(v1) == 0
                rem = rem | (rop(:, 1) == v1) | (rop(:, 2) == v1);
            end
            
            if hookConnDiff(v2) == 0
                rem = rem | (rop(:, 1) == v2) | (rop(:, 2) == v2);
            end

            rem = ~rem;

            rop = rop(rem, :);
            ropEdgeIds = ropEdgeIds(rem, :);
            edgeFilter = edgeFilter(rem, :);
        else
            rop = rop(2 : end, :);
            ropEdgeIds = ropEdgeIds(2 : end, :);
            edgeFilter = edgeFilter(2 : end, :);
        end
        
%         % TEST -> Just move the top lines to the bottom
%         rop = [rop(2 : end, :); rop(1, :)];
%         ropEdgeIds = [ropEdgeIds(2 : end, :); ropEdgeIds(1, :)];
%         edgeFilter = [edgeFilter(2 : end, :); edgeFilter(1, :)];
        
        %fprintf('hookConnDiff(51) == %d, inside == %d\n', hookConnDiff(51), any(rop(:) == 51));
    end
    
    assert(max(hookConnDiff) <= 1 && sum(hookConnDiff) <= 2);
    % fprintf('sum(hookConnDiff) == %d\n', sum(hookConnDiff));
    
    if sum(hookConnDiff) == 0
        startNode = usedNodeIds(randi(numUsedNodeIds));
    else
        ind = find(hookConnDiff);
        startNode = ind(1);
    end
    
%     if (size(rop, 1) == 1)
%         startNode = rop(1,1);
%     else
%         startNode = usedNodeIds(randi(numUsedNodeIds));
%     end
    
    fprintf('INFO: Added %d edges to graph to make it Eulerian\n', added);
    
    % Find Euler Path
    N = neighborMatrix(eCodes, x, aToBType, bToAType);
    
    INDEX = 1;
    FROM_DIR = 2;
    TO_DIR = 3;
    
    neighbors = [sum(N(:,:,1))' sum(N(:,:,2))' sum(N(:,:,3))' sum(N(:,:,4))'];

    if hookConnFromRight(startNode) > hookConnFromLeft(startNode)
        validConnTypes = [1 3]';
    elseif hookConnFromRight(startNode) < hookConnFromLeft(startNode)
        validConnTypes = [2 4]';
    else
        validConnTypes = [1 2 3 4]';
    end
    
    isCurrentSet = false;
    i = 1;
    while (~isCurrentSet)
        conn = validConnTypes(i);
        if neighbors(startNode, conn) ~= 0
            currentVertex = [startNode 0 conn];
            isCurrentSet = true;
        end
    end
    
    circuit = Stack(sum(x) + 1, 3);
    s = Stack(sum(x), 3);

    %rng(0);
    
    lastLeavingDir = zeros(numHooks, 1);
    
    while (neighbors(currentVertex(INDEX), currentVertex(TO_DIR)) > 0 || s.count() > 0)
        if (neighbors(currentVertex(INDEX), currentVertex(TO_DIR)) > 0)
            s.push(currentVertex);
            
            i = find(N(:, currentVertex(INDEX), currentVertex(TO_DIR)));
            
            neighborIndex = i(randi(size(i, 1)));
            
            % Update Neighbor Valence Information

            neighborFrom = currentVertex(TO_DIR);
            
            if currentVertex(TO_DIR) == 1
                neighborFrom = 2;
            elseif currentVertex(TO_DIR) == 2
                neighborFrom = 1;
            end
            
            if neighborFrom == 2 || neighborFrom == 4
                % Neighbor is reached from the LEFT -> Go to RIGHT (1 or 3)
                indices = [1 3]';
                %possibleEdges = indices;
                possibleEdges = indices(randperm(2));
                
                if neighbors(neighborIndex, possibleEdges(1)) > 0
                    neighborTo = possibleEdges(1);
                elseif neighbors(neighborIndex, possibleEdges(2)) > 0
                    neighborTo = possibleEdges(2);
                elseif lastLeavingDir(neighborIndex) > 0
                    neighborTo = lastLeavingDir(neighborIndex);
                else
                    neighborTo = possibleEdges(2);
                end
            else
                % Neighbor is reached from the RIGHT -> Go to LEFT (2 or 4)
                indices = [2 4]';
                %possibleEdges = indices;
                possibleEdges = indices(randperm(2));
                
                if neighbors(neighborIndex, possibleEdges(1)) > 0
                    neighborTo = possibleEdges(1);
                elseif neighbors(neighborIndex, possibleEdges(2)) > 0
                    neighborTo = possibleEdges(2);
                elseif lastLeavingDir(neighborIndex) > 0
                    neighborTo = lastLeavingDir(neighborIndex);
                else
                    neighborTo = possibleEdges(2);
                end
            end
            
            neighborVertex = [neighborIndex neighborFrom neighborTo];

            neighbors(currentVertex(INDEX), currentVertex(TO_DIR)) = ...
                neighbors(currentVertex(INDEX), currentVertex(TO_DIR)) - 1;
            
            neighbors(neighborVertex(INDEX), neighborVertex(FROM_DIR)) = ...
                neighbors(neighborVertex(INDEX), neighborVertex(FROM_DIR)) - 1;
            
            N(neighborVertex(INDEX), currentVertex(INDEX), currentVertex(TO_DIR)) = ...
                N(neighborVertex(INDEX), currentVertex(INDEX), currentVertex(TO_DIR)) - 1;
            
            N(currentVertex(INDEX), neighborVertex(INDEX), neighborVertex(FROM_DIR)) = ...
                N(currentVertex(INDEX), neighborVertex(INDEX), neighborVertex(FROM_DIR)) - 1;
            
            currentVertex = neighborVertex;
        else
            circuit.push(currentVertex);
            currentVertex = s.pop();
            
            lastLeavingDir(currentVertex(INDEX)) = currentVertex(TO_DIR);
            
            % Change dir if there is no neighbor
            if neighbors(currentVertex(INDEX), currentVertex(TO_DIR)) == 0
                if currentVertex(TO_DIR) == 1 && neighbors(currentVertex(INDEX), 3) > 0
                    currentVertex(TO_DIR) = 3;
                elseif currentVertex(TO_DIR) == 3 && neighbors(currentVertex(INDEX), 1) > 0
                    currentVertex(TO_DIR) = 1;
                elseif currentVertex(TO_DIR) == 2 && neighbors(currentVertex(INDEX), 4) > 0
                    currentVertex(TO_DIR) = 4;
                elseif currentVertex(TO_DIR) == 4 && neighbors(currentVertex(INDEX), 2) > 0
                    currentVertex(TO_DIR) = 2;
                end
            end
        end
    end
    
    circuit.push(currentVertex);
    
    data = flipud(circuit.data);
    
    stringList = zeros(size(data, 1) - 1, 4);
    stringList(:, 1) = data(1 : end - 1, 1);
    stringList(:, 2) = data(2 : end, 1);
    stringList(:, 3) = data(1 : end - 1, 3);
    stringList(:, 4) = data(2 : end, 2);
    
    hookSequence = data(:, 1);
    
    hookApproachLeftSide = data(2 : end, 2) == 2 | data(2 : end, 2) == 4;
    
    outsideMask = computeOutsideMask(stringList, auxNeigh);
end

function mask = computeOutsideMask(stringList, auxNeigh)
    mask = false(size(stringList, 1), 1);
    
    for i = 1 : size(stringList, 1)
        from = stringList(i, 1);
        to = stringList(i, 2);
        type = stringList(i, 3);
        if auxNeigh(to, from, type) > 0
            mask(i) = true;
            
            auxNeigh(to, from, type) = auxNeigh(to, from, type) - 1;
            
            if type == 3 || type == 4
                auxNeigh(from, to, type) = auxNeigh(from, to, type) - 1;
            elseif type == 1
                auxNeigh(from, to, 2) = auxNeigh(from, to, 2) - 1;
            elseif type == 2
                auxNeigh(from, to, 1) = auxNeigh(from, to, 1) - 1;
            end
        end
    end
end

function connected = isConnected(eCodes, x)
    % Map Hook Numbers to virtual extended hook number range
    numHooks = max(max(eCodes));
    numEdges = size(eCodes, 1);
    numHookConn = 4;
    
    numPointToPointConnections = numEdges / numHookConn;
    
    extOff = numHooks * repmat((0 : 1 : (numHookConn - 1))', numPointToPointConnections, 1);
    
    extEdgeCodes = eCodes + [extOff extOff];
    extEdgeCodes = extEdgeCodes(x == 1, :);
    
    % Add edges to connect all virtual sub-nodes in each hook
    %numExtEdges = numEdges + numHooks * numHookConn;
    
    % (1+0*256) (1+1*256)
    % (1+1*256) (1+2*256)
    % (1+2*256) (1+3*256)
    
    % (2+0*256) (2+1*256)
    % (2+1*256) (2+2*256)
    % (2+2*256) (2+3*256)
    base = repelem((1 : 1 : numHooks)', numHookConn - 1);
    add = 256 * repmat((0 : 1 : 2)', numHooks, 1);
    addEdges = [(base + add) (base + add + numHooks)];
    extEdgeCodes = [extEdgeCodes; addEdges];
    
    tempG = graph(extEdgeCodes(:, 1), extEdgeCodes(:, 2));
    bins = conncomp(tempG);
    
    uniqueBins = unique(bins);
    
    connected = (numel(uniqueBins) == 1);
end

function N = neighborMatrix(eCodes, x, aToBType, bToAType)
    numHooks = max(max(eCodes));
    data = [[eCodes; fliplr(eCodes)] [x; x] [aToBType; bToAType]];
    data = data(data(:, 3) > 0, :);
    
    mask1 = data(:, 4) == 1;
    mask2 = data(:, 4) == 2;
    mask3 = data(:, 4) == 3;
    mask4 = data(:, 4) == 4;
    
    N = cat(3, ...
        full(sparse(data(mask1, 2), data(mask1, 1), data(mask1, 3), numHooks, numHooks)), ...
        full(sparse(data(mask2, 2), data(mask2, 1), data(mask2, 3), numHooks, numHooks)), ...
        full(sparse(data(mask3, 2), data(mask3, 1), data(mask3, 3), numHooks, numHooks)), ...
        full(sparse(data(mask4, 2), data(mask4, 1), data(mask4, 3), numHooks, numHooks)));
end