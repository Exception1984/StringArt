classdef GreedyMultiSamplingDataObjectL2 < handle

    properties       
        AEdgeIndicesToPixelCodes
        
        highResEdgePixelIndices
        highResEdgePixelValues
        highResCorrespEdgeIndices
        highResIndexToIndexMap
        
        lowResEdgePixelIndices
        lowResEdgePixelValues
        lowResCorrespEdgeIndices
        lowResIndexToIndexMap
        
        highResReIndexMap
        lowResReIndexMap
        
        highResTransposeIndexMap
        lowResTransposeIndexMap
        
        filterWeight
        highResToLowResMap
        lowResToHighResMap
        
        matrixPath
        
        corrMap
        bNativeRes
        EdgeToHook
        HookToEdge
        
        diffToBlankSquaredErrors
        diffToBlankSquaredErrorSum
        
        highResImageWidth
        lowResImageWidth

        importanceWeight
        
        currentRecon
        currentReconUnclamped
        currentReconNativeRes
        
        currentReconSquare
        
        numEdges
        numHooks
        
        allIndices
        
        reachablePixelsMask
        reachablePixelsMaskNativeRes
        
        x
        hookCount
        pickedEdgesSequence
        
        FAdding
        FRemoving
        rmseValue
        
        consecutive
        incident
        
        stringList
        numStrings
        
        minCircleLength
        maxNumEdgeUsage
        maxNumHookUsage
        maxNumStrings
        
        latelyVisitedPins
        
        currentHook
        
        validEdgesMask
        
        removalMode
        removableEdgeIndices
        
        numLeftEdgesPerHook
        numRightEdgesPerHook
        hookToLeftEdges
        hookToRightEdges
        hookSideBalance
        
        importanceMap
        
        clamping
    end
    
    methods
        function obj = GreedyMultiSamplingDataObjectL2(img, domainWidth, windowSize, maxNumStrings, minAngle, numHooks, importanceMap, matrixPath)
            [obj.EdgeToHook, obj.HookToEdge] = edgeCodes(numHooks);
            
            obj.clamping = true;
            
            obj.hookSideBalance = true;
            
            obj.matrixPath = matrixPath;
            
            obj.createImageVector(img, domainWidth, windowSize);
            
            if nargin < 8 || isempty(importanceMap)
                importanceMap = ones(domainWidth, domainWidth);
            end
            
            obj.createImportanceMapVector(importanceMap, domainWidth);
            
            obj.lowResImageWidth = domainWidth;
            obj.highResImageWidth = domainWidth * windowSize;
            
            obj.computeTransposeIndexMaps();
            
            obj.currentRecon = zeros(obj.highResImageWidth * obj.highResImageWidth, 1);
            obj.currentReconUnclamped = obj.currentRecon;
            obj.currentReconSquare = zeros(obj.highResImageWidth, obj.highResImageWidth);
            
            obj.currentReconNativeRes = zeros(obj.lowResImageWidth * obj.lowResImageWidth, 1);
            
            obj.numEdges = size(obj.EdgeToHook, 1);
            obj.numHooks = numHooks;
            
            obj.validEdgesMask = true(obj.numEdges, 1);
            
            obj.validEdgesMask = loadFabricatableEdgesMask(obj.matrixPath);
            assert(obj.numEdges == numel(obj.validEdgesMask));

            obj.numLeftEdgesPerHook = zeros(obj.numHooks, 1);
            obj.numRightEdgesPerHook = zeros(obj.numHooks, 1);

            obj.initLeftAndRighEdges();
            
            if minAngle > 0
                obj.computeValidEdgesMask(minAngle);
            end
            
            obj.loadAIndexMatrices();
            
            obj.reachablePixelsMask = false((obj.highResImageWidth * obj.highResImageWidth), 1);
            obj.reachablePixelsMask(obj.highResEdgePixelIndices) = true;
            
            obj.reachablePixelsMaskNativeRes = false(obj.lowResImageWidth * obj.lowResImageWidth, 1);
            obj.reachablePixelsMaskNativeRes(obj.lowResEdgePixelIndices) = true;
            
            if minAngle > 0
                obj.computeCoreZonePixels(minAngle);
            end
            
            obj.filterWeight = 1.0 / (windowSize * windowSize);
            
            obj.highResReIndexMap = zeros(size(obj.corrMap, 2), 1);
            obj.lowResReIndexMap = zeros(size(obj.corrMap, 1), 1);

            obj.removalMode = false;
            
            obj.allIndices = (1 : size(obj.corrMap, 2))';
            
            obj.x = zeros(obj.numEdges, 1);
            obj.removableEdgeIndices = zeros(0, 1);
            obj.hookCount = zeros(size(obj.HookToEdge,1), 1);
            obj.pickedEdgesSequence = zeros(0, 1);
            
            obj.consecutive = false;
            obj.incident = find(obj.validEdgesMask);
            
            obj.maxNumStrings = maxNumStrings;
            obj.stringList = zeros(0, 3);
            obj.numStrings = 1;
            
            obj.minCircleLength = 1;
            obj.maxNumEdgeUsage = 1;
            obj.maxNumHookUsage = 100000;
            
            obj.currentHook = 0;
            
            obj.initStateVectors();
            obj.initLatelyVisitedPins();
            obj.initIndexMaps();
        end
        
        function initLeftAndRighEdges(obj)
            obj.hookToLeftEdges = zeros(obj.numHooks, 2 * (obj.numHooks - 1));
            obj.hookToRightEdges = zeros(obj.numHooks, 2 * (obj.numHooks - 1));
            
            for i = 1 : obj.numHooks
                data = obj.HookToEdge{i};
                leftMask = data(:, 2) == 2 | data(:, 2) == 4;
                rightMask = ~leftMask;
                
                obj.hookToLeftEdges(i, :) = data(leftMask, 1)';
                obj.hookToRightEdges(i, :) = data(rightMask, 1)';
            end
        end
        
        function computeTransposeIndexMaps(obj)
            highResIndices = reshape((1 : obj.highResImageWidth * obj.highResImageWidth)', [obj.highResImageWidth, obj.highResImageWidth])';
            obj.highResTransposeIndexMap = highResIndices(:);
            
            lowResIndices = reshape((1 : obj.lowResImageWidth * obj.lowResImageWidth)', [obj.lowResImageWidth , obj.lowResImageWidth])';
            obj.lowResTransposeIndexMap = lowResIndices(:);
        end
        
        function createImageVector(obj, img, domainWidth, windowSize)
            if size(img, 1) ~= domainWidth || size(img, 2) ~= domainWidth
                img = imresize(img, [domainWidth, domainWidth]);
            end
            
            obj.bNativeRes = full(buildImageVector(img, domainWidth));
            obj.corrMap = multiSampleCorrespondenceMap(domainWidth, windowSize);
        end
        
        function createImportanceMapVector(obj, importanceMap, domainWidth)
            if size(importanceMap, 1) ~= domainWidth || size(importanceMap, 2) ~= domainWidth
                importanceMap = imresize(importanceMap, [domainWidth, domainWidth]);
            end
            
            obj.importanceMap = full(buildImageVector(importanceMap, domainWidth));
        end
        
        function initIndexMaps(obj)
            if (obj.lowResImageWidth == obj.highResImageWidth)
                obj.highResToLowResMap = uint32((1 : obj.highResImageWidth * obj.highResImageWidth)');
                obj.lowResToHighResMap = uint32(obj.highResToLowResMap);
            else
                obj.highResToLowResMap = uint32(highResToLowResIndex((1 : obj.highResImageWidth * obj.highResImageWidth)', obj.highResImageWidth, obj.lowResImageWidth));
                obj.lowResToHighResMap = uint32(lowResToHighResIndex((1 : obj.lowResImageWidth * obj.lowResImageWidth)', obj.lowResImageWidth, obj.highResImageWidth));
            end
        end
        
        function loadAIndexMatrices(obj)
            validEdgeIndices = find(obj.validEdgesMask);
            invalidEdgeIndices = find(~obj.validEdgesMask);
            
            % obj.lowResImageWidth = domainWidth;
            % obj.highResImageWidth = domainWidth * windowSize;
            
            % High Resolution Version
            tic
            fprintf('Loading High Resolution Index Matrices... ');
            [obj.highResEdgePixelIndices, obj.highResEdgePixelValues, obj.highResCorrespEdgeIndices, obj.highResIndexToIndexMap] = loadAIndexMatrices(obj.matrixPath, obj.highResImageWidth, invalidEdgeIndices);
            fprintf(['DONE (' num2str(toc) ' s)\n']);
            
            % Low Resolution Version
            tic
            fprintf('Loading Low Resolution Index Matrices... ');
            [obj.lowResEdgePixelIndices, obj.lowResEdgePixelValues, obj.lowResCorrespEdgeIndices, obj.lowResIndexToIndexMap] = loadAIndexMatrices(obj.matrixPath, obj.lowResImageWidth, invalidEdgeIndices);
            fprintf(['DONE (' num2str(toc) ' s)\n']);
            
            tic
            fprintf('Computing Edge-Index To Pixel-Index Matrix... ');
            obj.AEdgeIndicesToPixelCodes = cell(obj.numEdges, 2);
            
            tempMatrix = sparse(double(obj.highResEdgePixelIndices), double(obj.highResCorrespEdgeIndices), obj.highResEdgePixelValues, obj.highResImageWidth * obj.highResImageWidth, obj.numEdges);
            
            for j = 1 : numel(validEdgeIndices)
                i = validEdgeIndices(j);
                [row, ~, val] = find(tempMatrix(:, i));
                obj.AEdgeIndicesToPixelCodes{i, 1} = uint32(row);
                obj.AEdgeIndicesToPixelCodes{i, 2} = val;
            end
            
            fprintf(['DONE (' num2str(toc) ' s)\n']);
            
            clear tempMatrix;
        end
        
        function computeValidEdgesMask(obj, minAngle)
            %hooks = getUsedHookPositions();
            
            hooks = matrixFromFile(fullfile(obj.matrixPath, 'pinPositions'), 'double');
    
            r = sqrt(hooks(:,1).* hooks(:,1) + hooks(:,2) .* hooks(:,2));
            hooks = hooks ./ [r r];

            % Filter out edges that connect very close hooks
            hooks_x = hooks(:,1);
            hooks_y = hooks(:,2);

            x1 = hooks_x(obj.EdgeToHook(:,1));
            y1 = hooks_y(obj.EdgeToHook(:,1));
            x2 = hooks_x(obj.EdgeToHook(:,2));
            y2 = hooks_y(obj.EdgeToHook(:,2));

            edgeAngles = acos(x1 .* x2 + y1 .* y2);

            obj.validEdgesMask(edgeAngles < minAngle) = false;
        end
        
        function fac = computeCoreZonePixels(obj, minAngle)
            fac = coreDomainFactor(minAngle);
            
            if fac < 1
                w = obj.lowResImageWidth * obj.lowResImageWidth;
                a = sqrt(w);
                mask = true(a, a);
                mask = maskImage(mask, 0.5 * fac * a);
                
                obj.reachablePixelsMaskNativeRes = and(mask(:), obj.reachablePixelsMaskNativeRes);
            end
        end
        
        function initStateVectors(obj)
            obj.diffToBlankSquaredErrors = (obj.importanceMap .* obj.bNativeRes).^2;
            obj.diffToBlankSquaredErrorSum = sum(obj.diffToBlankSquaredErrors);

            obj.rmseValue = sqrt(obj.diffToBlankSquaredErrorSum / numel(obj.bNativeRes));
            
            importancePerEdgePixel = obj.importanceMap(obj.lowResEdgePixelIndices);
            
            sumOfSquaredErrorsPerEdgeAdding = accumarray(obj.lowResCorrespEdgeIndices, ...
                (importancePerEdgePixel .* (obj.bNativeRes(obj.lowResEdgePixelIndices) - obj.lowResEdgePixelValues)).^2, [obj.numEdges, 1]);
            
            sumOfSquaredErrorsPerEdgeRemoving = accumarray(obj.lowResCorrespEdgeIndices, ...
                (importancePerEdgePixel .* (obj.bNativeRes(obj.lowResEdgePixelIndices) + obj.lowResEdgePixelValues)).^2, [obj.numEdges, 1]);
            
            diffToBlankSumPerEdge = accumarray(obj.lowResCorrespEdgeIndices, obj.diffToBlankSquaredErrors(obj.lowResEdgePixelIndices), [obj.numEdges, 1]);
            
            obj.FAdding = obj.diffToBlankSquaredErrorSum - diffToBlankSumPerEdge + sumOfSquaredErrorsPerEdgeAdding;
            obj.FRemoving = obj.diffToBlankSquaredErrorSum - diffToBlankSumPerEdge + sumOfSquaredErrorsPerEdgeRemoving;
        end
        
        function initLatelyVisitedPins(obj)
            obj.latelyVisitedPins = zeros(1, 0);
        end

        function [val, i] = findBestString(obj)
            if obj.removalMode
                [val, j] = min(obj.FRemoving(obj.removableEdgeIndices));
                %ind = find(obj.removableEdgeIndices);
                 i = obj.removableEdgeIndices(j); %ind(j);
            else
                [val, j] = min(obj.FAdding(obj.incident));
                %ind = find(obj.incident);
                i = obj.incident(j); %ind(j);
            end
            
            fprintf('\tF1 when picking edge Nr. %i: %16.16f\n', i, val);
        end
        
        function val = F1(obj, i)
            val = obj.FRemoving(i);
            %fprintf('F1 when picking edge Nr. %i: %16.16f\n', i, val);
        end
        
        function val = F2(obj)            
            val = sum((obj.bNativeRes - obj.currentReconNativeRes).^2);
        end
        
        function chooseStringAndUpdate(obj, i)
            % Find all relevant indices
            edgePixelIndices = obj.AEdgeIndicesToPixelCodes{i, 1};
            edgeValues = obj.AEdgeIndicesToPixelCodes{i, 2};
            nativeResIndices = unique(obj.highResToLowResMap(edgePixelIndices));

            highResIndices = obj.lowResToHighResMap(nativeResIndices, :);
            highResIndices = highResIndices(:);
            
            %preUpdateHighResRecon = obj.currentRecon(highResIndices);
            preUpdateHighResReconUnclamped = obj.currentReconUnclamped(highResIndices);
            preUpdateLowResRecon = obj.currentReconNativeRes(nativeResIndices);
            
            dif = 1;
            
            if obj.removalMode
                mask = obj.removableEdgeIndices == i;
                if ~any(mask)
                    error('Edge %i can not be removed.\n', i);
                else
                    dif = -1;
                    
                    if (obj.x(i) - 1 == 0)
                        obj.removableEdgeIndices = obj.removableEdgeIndices(~mask);
                    end
                    
                    mask = obj.stringList(:, 1) == i;
                    ind = find(mask);
                    mask(:) = true;
                    mask(ind(1)) = false;
                    obj.stringList = obj.stringList(mask, :);
                    obj.pickedEdgesSequence = obj.pickedEdgesSequence(mask, :);
                end
            else
                obj.removableEdgeIndices = [obj.removableEdgeIndices; i];
                obj.stringList = [obj.stringList; i 0 0]; %(obj.numStrings, 1) = i;
                obj.pickedEdgesSequence = [obj.pickedEdgesSequence' i]';
            end
            
            % Update data structures
            obj.x(i) = obj.x(i) + dif;

            % Update Current Recon
            obj.currentReconUnclamped(edgePixelIndices) = obj.currentReconUnclamped(edgePixelIndices) + dif * edgeValues;
            
            if (obj.clamping)
                obj.currentRecon(edgePixelIndices) = min(1.0, obj.currentReconUnclamped(edgePixelIndices));
            else
                obj.currentRecon(edgePixelIndices) = obj.currentReconUnclamped(edgePixelIndices);
            end
            

            %obj.currentReconSquare(obj.highResTransposeIndexMap(edgePixelIndices)) = obj.currentRecon(edgePixelIndices);
            %obj.showCurrent();
            %imshow(obj.currentReconSquare);

            obj.currentReconNativeRes(nativeResIndices) = obj.corrMap(nativeResIndices, :) * obj.currentRecon;

            % testCurrentReconNativeRes = imresize(obj.currentReconSquare, [obj.lowResImageWidth, obj.lowResImageWidth], 'box');

            fprintf('\tF2 when picking edge Nr. %i: %16.16f\n\n', i, sum((obj.importanceMap .* (obj.bNativeRes - obj.currentReconNativeRes)).^2));
            %obj.showCurrent();

            preUpdateErrors = obj.diffToBlankSquaredErrors(nativeResIndices);
            obj.diffToBlankSquaredErrorSum = obj.diffToBlankSquaredErrorSum - sum(preUpdateErrors);
            obj.diffToBlankSquaredErrors(nativeResIndices) = (obj.importanceMap(nativeResIndices) .* (obj.bNativeRes(nativeResIndices) - obj.currentReconNativeRes(nativeResIndices))).^2;
            postUpdateErrors = obj.diffToBlankSquaredErrors(nativeResIndices);
            obj.diffToBlankSquaredErrorSum = obj.diffToBlankSquaredErrorSum + sum(postUpdateErrors);
            obj.rmseValue = sqrt(obj.diffToBlankSquaredErrorSum / numel(obj.bNativeRes));

            obj.updateEdgeErrors(nativeResIndices, highResIndices, preUpdateLowResRecon, preUpdateHighResReconUnclamped, preUpdateErrors, postUpdateErrors);

            % Update incidence vector
            obj.updateIncidenceVector(i);

            % Increment / Decrement Num Strings
            obj.numStrings = obj.numStrings + dif;
        end
        
        function updateEdgeErrors(obj, lowResIndices, highResIndices, preUpdateLowResRecon, preUpdateHighResReconUnclamped, preUpdateErrors, postUpdateErrors)
            % Update non-intersecting pixel positions
            pre  = sum(preUpdateErrors);
            post = sum(postUpdateErrors);
            
            % First, falsly update all edges and fix intersection errors afterwards
            obj.FAdding(:) = obj.FAdding(:) - pre + post;
            obj.FRemoving(:) = obj.FRemoving(:) - pre + post;
            
            secMask = any(obj.lowResIndexToIndexMap(:, lowResIndices), 2);
            truncHighResIndicesMask = highResIndices <= size(obj.highResIndexToIndexMap, 2);
            truncHighResIndices = highResIndices(truncHighResIndicesMask);
            highResSecMask = any(obj.highResIndexToIndexMap(:, truncHighResIndices), 2);
            
            secEdgePixInd = obj.lowResEdgePixelIndices(secMask);
            
            obj.lowResReIndexMap(lowResIndices) = (1 : numel(lowResIndices))';
            reIndices = obj.lowResReIndexMap(secEdgePixInd);
            
            preAtIndices = preUpdateErrors(reIndices);
            postAtIndices = postUpdateErrors(reIndices);
            
            secCorrEdgeInd = obj.lowResCorrespEdgeIndices(secMask);
            
            preCorr = accumarray(secCorrEdgeInd, preAtIndices, [obj.numEdges 1]);
            postCorr = accumarray(secCorrEdgeInd, postAtIndices,  [obj.numEdges 1]);
            
            % Fix intersection errors
            obj.FAdding(:) = obj.FAdding(:) - postCorr + preCorr;
            obj.FRemoving(:) = obj.FRemoving(:) - postCorr + preCorr;
            
            % Update intersecting pixel positions
            highResSecEdgePixInd = obj.highResEdgePixelIndices(highResSecMask);
            highResSecEdgePixVal = obj.highResEdgePixelValues(highResSecMask);
            highResCorrEdgeInd = obj.highResCorrespEdgeIndices(highResSecMask);
            highResToLowResInd = obj.highResToLowResMap(highResSecEdgePixInd);
            
            m = double(max(highResToLowResInd));
            
            outputIDs = rowAndColToIndex(double(highResCorrEdgeInd), double(highResToLowResInd), m);
            uniqueOutputIDs = unique(outputIDs);
            
            [filteredCorrEdgeInd, filteredPixInd] = indexToRowAndCol(uniqueOutputIDs, m);
            filteredReIndices = obj.lowResReIndexMap(filteredPixInd);
            
            numIDs = numel(uniqueOutputIDs);
            outputReIndex = (1 : numIDs)';
            [Lia,Locb] = ismember(outputIDs, uniqueOutputIDs);
            outputReIndex = outputReIndex(Locb);
            
            obj.highResReIndexMap(highResIndices) = (1 : numel(highResIndices))';
            highResReIndices = obj.highResReIndexMap(highResSecEdgePixInd);
            
            %% PRE UPDATE
            if obj.clamping
                preUpdateHighResRecon = min(1.0, preUpdateHighResReconUnclamped);
            else
                preUpdateHighResRecon = preUpdateHighResReconUnclamped;
            end

            preUpdateHighResVal = preUpdateHighResRecon(highResReIndices);

            preUpdateHighResValUnclamped = preUpdateHighResReconUnclamped(highResReIndices);

            if obj.clamping
                preUpdateHighResValMinusEdges = min(1.0, preUpdateHighResValUnclamped - highResSecEdgePixVal);
                preUpdateHighResValPlusEdges = min(1.0, preUpdateHighResValUnclamped + highResSecEdgePixVal);
            else
                preUpdateHighResValMinusEdges = preUpdateHighResValUnclamped - highResSecEdgePixVal;
                preUpdateHighResValPlusEdges = preUpdateHighResValUnclamped + highResSecEdgePixVal;
            end

            % Filter pre update highRes edges
            filteredPreUpdateHighResVal = obj.filterWeight * accumarray(outputReIndex, preUpdateHighResVal);
            filteredPreUpdateHighResValMinusEdges = obj.filterWeight * accumarray(outputReIndex, preUpdateHighResValMinusEdges);
            filteredPreUpdateHighResValPlusEdges = obj.filterWeight * accumarray(outputReIndex, preUpdateHighResValPlusEdges);
            
            failurePreUpdatePerEdgeAdding = accumarray(filteredCorrEdgeInd, ...
                (obj.importanceMap(filteredPixInd) .* (obj.bNativeRes(filteredPixInd) - (preUpdateLowResRecon(filteredReIndices) - filteredPreUpdateHighResVal + filteredPreUpdateHighResValPlusEdges))).^2,  [obj.numEdges 1]);
            
            failurePreUpdatePerEdgeRemoving = accumarray(filteredCorrEdgeInd, ...
                (obj.importanceMap(filteredPixInd) .* (obj.bNativeRes(filteredPixInd) - (preUpdateLowResRecon(filteredReIndices) - filteredPreUpdateHighResVal + filteredPreUpdateHighResValMinusEdges))).^2,  [obj.numEdges 1]);
            
            %% POST UPDATE
            postUpdateHighResRecon = obj.currentRecon(highResIndices);
            postUpdateHighResVal = postUpdateHighResRecon(highResReIndices);

            postUpdateHighResReconUnclamped = obj.currentReconUnclamped(highResIndices);
            postUpdateHighResValUnclamped = postUpdateHighResReconUnclamped(highResReIndices);
            
            if obj.clamping
                postUpdateHighResValMinusEdges = min(1.0, postUpdateHighResValUnclamped - highResSecEdgePixVal);
                postUpdateHighResValPlusEdges = min(1.0, postUpdateHighResValUnclamped + highResSecEdgePixVal);
            else
                postUpdateHighResValMinusEdges = postUpdateHighResValUnclamped - highResSecEdgePixVal;
                postUpdateHighResValPlusEdges = postUpdateHighResValUnclamped + highResSecEdgePixVal;
            end

            postUpdateLowResRecon = obj.currentReconNativeRes(lowResIndices);
            
            % Filter post update highRes edges
            filteredPostUpdateHighResVal = obj.filterWeight * accumarray(outputReIndex, postUpdateHighResVal);
            filteredPostUpdateHighResValMinusEdges = obj.filterWeight * accumarray(outputReIndex, postUpdateHighResValMinusEdges);
            filteredPostUpdateHighResValPlusEdges = obj.filterWeight * accumarray(outputReIndex, postUpdateHighResValPlusEdges);
            
            failurePostUpdatePerEdgeAdding = accumarray(filteredCorrEdgeInd, ...
                (obj.importanceMap(filteredPixInd) .* (obj.bNativeRes(filteredPixInd) - (postUpdateLowResRecon(filteredReIndices) - filteredPostUpdateHighResVal + filteredPostUpdateHighResValPlusEdges))).^2,  [obj.numEdges 1]);
            
            failurePostUpdatePerEdgeRemoving = accumarray(filteredCorrEdgeInd, ...
                (obj.importanceMap(filteredPixInd) .* (obj.bNativeRes(filteredPixInd) - (postUpdateLowResRecon(filteredReIndices) - filteredPostUpdateHighResVal + filteredPostUpdateHighResValMinusEdges))).^2,  [obj.numEdges 1]);
            
            obj.FAdding(:) = obj.FAdding(:) - failurePreUpdatePerEdgeAdding + failurePostUpdatePerEdgeAdding;
            obj.FRemoving(:) = obj.FRemoving(:) - failurePreUpdatePerEdgeRemoving + failurePostUpdatePerEdgeRemoving;
        end
        
        function updateIncidenceVector(obj, edgeID)
            dif = 1;
            
            if obj.removalMode
                dif = -1;
            end
            
            edgeType = mod(edgeID, 4);
            
            if edgeType == 0
                edgeType = 4;
            end
            
            hooks = obj.EdgeToHook(edgeID, :)';
            obj.hookCount(hooks(1)) = obj.hookCount(hooks(1)) + dif;
            obj.hookCount(hooks(2)) = obj.hookCount(hooks(2)) + dif;
            
            hookA = hooks(1);
            hookB = hooks(2);
            
            if edgeType == 1
                % RIGHT for hookA, LEFT for hookB
                obj.numRightEdgesPerHook(hookA) = obj.numRightEdgesPerHook(hookA) + dif;
                obj.numLeftEdgesPerHook(hookB) = obj.numLeftEdgesPerHook(hookB) + dif;
            elseif edgeType == 2
                % LEFT for hookA, RIGHT for hookB
                obj.numLeftEdgesPerHook(hookA) = obj.numLeftEdgesPerHook(hookA) + dif;
                obj.numRightEdgesPerHook(hookB) = obj.numRightEdgesPerHook(hookB) + dif;
            elseif edgeType == 3
                % RIGHT for hookA, RIGHT for hookB
                obj.numRightEdgesPerHook(hookA) = obj.numRightEdgesPerHook(hookA) + dif;
                obj.numRightEdgesPerHook(hookB) = obj.numRightEdgesPerHook(hookB) + dif;
            else % edgeType == 4
                % LEFT for hookA, LEFT for hookB
                obj.numLeftEdgesPerHook(hookA) = obj.numLeftEdgesPerHook(hookA) + dif;
                obj.numLeftEdgesPerHook(hookB) = obj.numLeftEdgesPerHook(hookB) + dif;
            end
            
            if (obj.consecutive == true)
                if (obj.currentHook == 0)
                    
                    % This was the first edge
                    edgesToH1 = obj.HookToEdge{hooks(1)};
                    edgesToH2 = obj.HookToEdge{hooks(2)};
                    obj.incident = [edgesToH1; edgesToH2];
                    
                    [val, i] = obj.findBestString();
                
                    if (i <= size(edgesToH1, 1))
                        obj.currentHook = hooks(2);
                    else
                        obj.currentHook = hooks(1);
                    end
                    
                    obj.latelyVisitedPins = [obj.latelyVisitedPins' obj.currentHook]';
                end
                
                % Not the first edge
                obj.stringList(obj.numStrings, 2) = obj.currentHook;

                if (hooks(1) == obj.currentHook)
                    obj.currentHook = hooks(2);
                else
                    obj.currentHook = hooks(1);
                end

                obj.stringList(obj.numStrings, 3) = obj.currentHook;
                obj.hookCount(obj.currentHook) = obj.hookCount(obj.currentHook) + 1;

                obj.latelyVisitedPins = [obj.latelyVisitedPins' obj.currentHook]';

                if size(obj.latelyVisitedPins, 1) == obj.minCircleLength
                    obj.latelyVisitedPins = obj.latelyVisitedPins(2 : end);
                end
                
                obj.incident = obj.HookToEdge{obj.currentHook};
                
                % Remove incident edges that would create a circle < minCircleLength
                [maxUsedHooks, v] = find(obj.hookCount == obj.maxNumHookUsage);
                illegalPins = union(obj.latelyVisitedPins, maxUsedHooks);
                obj.incident = setdiff(obj.incident, obj.computeIllegalEdgeIndices(obj.currentHook, illegalPins));
            else
                % Non-consecutive string path
                
                % Update Incidence
                % Newly construct the incidence vector
                obj.incident = obj.validEdgesMask;
                
                % Remove overused edges
                obj.incident(obj.x >= obj.maxNumEdgeUsage) = false;
                
                % Remove edges that belong to overused hooks
                overusedHooksEdgeIndices = cell2mat(obj.HookToEdge(obj.hookCount > obj.maxNumHookUsage, 1));
                
                if ~isempty(overusedHooksEdgeIndices)
                    obj.incident(overusedHooksEdgeIndices(:, 1)) = false;
                end
                
                % Remove edges that would lead to an inbalance of the
                % sides of the hooks
                if obj.hookSideBalance
                    removeLeftSideEdgesMask  = obj.numLeftEdgesPerHook - obj.numRightEdgesPerHook > 0;
                    removeRightSideEdgesMask = obj.numRightEdgesPerHook - obj.numLeftEdgesPerHook > 0;

                    obj.incident(obj.hookToLeftEdges(removeLeftSideEdgesMask, :)) = false;
                    obj.incident(obj.hookToRightEdges(removeRightSideEdgesMask, :)) = false;
                end

                obj.incident = find(obj.incident);
                
                % Update removable edge indices
                obj.removableEdgeIndices = obj.x > 0;
                
                if obj.hookSideBalance
                    removeLeftSideEdgesMask = obj.numLeftEdgesPerHook - obj.numRightEdgesPerHook < 0;
                    removeRightSideEdgesMask = obj.numRightEdgesPerHook - obj.numLeftEdgesPerHook < 0;
                    obj.removableEdgeIndices (obj.hookToLeftEdges(removeLeftSideEdgesMask, :)) = false;
                    obj.removableEdgeIndices (obj.hookToRightEdges(removeRightSideEdgesMask, :)) = false;
                end
                
                obj.removableEdgeIndices = find(obj.removableEdgeIndices);
                
%                 % Remove incident edges that would overuse a hook
%                 if ~obj.removalMode
%                     [maxUsedHooks, v] = find(obj.hookCount == obj.maxNumHookUsage);
%                     obj.incident = setdiff(obj.incident, obj.computeIllegalEdgeIndices(hooks(1), maxUsedHooks));
%                     obj.incident = setdiff(obj.incident, obj.computeIllegalEdgeIndices(hooks(2), maxUsedHooks));
%                 end
            end
            
%             if obj.removalMode
%                 obj.incident = union(obj.incident, edgeID);
%             else
%                 %Remove overused edges
%                 obj.removeOverusedEdgesFromIncidenceVector();
%             end
        end
        
        function indices = computeIllegalEdgeIndices(obj, hook, illegalPins)
            if (hook == illegalPins)
                indices = zeros(1,0);
                return;
            end
            
            if (size(illegalPins, 1) == 1)
                illegalPins = [illegalPins illegalPins]';
            end

            lateleyVisitedIndices = repmat(illegalPins', size(obj.EdgeToHook, 1), 1);
            latelyFrom = any((lateleyVisitedIndices == repmat(obj.EdgeToHook(:, 1), 1, size(illegalPins, 1)))')';
            latelyTo = any((lateleyVisitedIndices == repmat(obj.EdgeToHook(:, 2), 1, size(illegalPins, 1)))')';

            curr = repmat(hook, size(obj.EdgeToHook, 1), 1);
            currFrom = curr == repmat(obj.EdgeToHook(:, 1), 1, 1);
            currTo = curr == repmat(obj.EdgeToHook(:, 2), 1, 1);

            res = or(and(latelyFrom, currTo), and(latelyTo, currFrom));
            [k, l, ~] = find(res);

            indices = k;
        end
        
        function removeOverusedEdgesFromIncidenceVector(obj)
            [maxUsedEdges, v] = find(obj.x >= obj.maxNumEdgeUsage);
            obj.incident = setdiff(obj.incident, maxUsedEdges);
        end
        
        function nccValue = getNccValue(obj)
            nccValue = imageCorrelation(obj.currentReconNativeRes, obj.bNativeRes);
        end
        
        function setImportanceMap (obj, importanceMap)
            im = importanceMap(obj.reachablePixelsMaskNativeRes);
            obj.importanceWeight(:) = im(:);
        end
        
        function showCurrent(obj)
            figure(1);
            imshow(imresize(1 - flipud(reshape(obj.currentRecon, [obj.highResImageWidth, obj.highResImageWidth])'), [512, 512]));
        end
        
        function rmse = getRmseValue(obj)
            rmse = obj.rmseValue;
        end
        
        function setRemovalMode(obj, mode)
            if mode ~= obj.removalMode
                obj.removalMode = mode;

                if obj.removalMode
                    if obj.consecutive
                        error('Removing Edges in a consecutive string chain is not allowed!');
                    end
                end
            end            
        end
    end
end