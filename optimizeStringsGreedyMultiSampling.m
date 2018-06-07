function [x, pickedEdgesSequence] = optimizeStringsGreedyMultiSampling(img, domainWidth, windowSize, minAngle, numHooks, importanceMap, matrixPath)
    obj = GreedyMultiSamplingDataObjectL2(img, domainWidth, windowSize, 100000, minAngle, numHooks, importanceMap, matrixPath);
    
    obj.consecutive = false;
    obj.maxNumEdgeUsage = 1;
    
    iterativeStepSize = 1;
    iteration = 0;
    bestRmseValue = inf;
    bestNumEdges = 0;
    numBadRuns = 0;

    while(iteration < obj.maxNumStrings)
        iteration = iteration + 1;
        
        fprintf('Iteration Nr. %d\n', iteration);
        
        [m, i] = obj.findBestString();
        
        if isempty(i)
            break;
        end
        
        if (size(i, 1) > 1)
            i = datasample(i,1);
        end
        
        obj.chooseStringAndUpdate(i);
        
        rmse = obj.getRmseValue();

        if (rmse < bestRmseValue)
            bestRmseValue = rmse;
            bestNumEdges = iteration;
            numBadRuns = 0;
        else
            numBadRuns = numBadRuns + 1;
        end
        
        if numBadRuns == 1000
            break;
        end
    end

    pureGreedyBestRmse = bestRmseValue;
    pickedEdgesSequence = obj.pickedEdgesSequence(1:bestNumEdges);
    
    if iterativeStepSize > 0
        % Try to improve result by iterative greedy approach
        
        %Remove 'overshoot'
        removeOvershoot(obj, numel(obj.pickedEdgesSequence) - bestNumEdges);
        
        %Iterative greedy approach
        numBadRuns = 0;
        modes = [true false]';
        
        removedEdgeIndices = zeros(iterativeStepSize, 1);
        addedEdgeIndices = zeros(iterativeStepSize, 1);
        
        while numBadRuns < 1000
            for k = 1 : 2
                if modes(k)
                    fprintf('Iterative Removal...\n');
                else
                    fprintf('Iterative Addition...\n');
                end
                
                obj.setRemovalMode(modes(k));
                for s = 1 : iterativeStepSize
                    [m, i] = obj.findBestString();
                    
                    if isempty(i)
                        break;
                    end

                    if (size(i, 1) > 1)
                        i = datasample(i,1);
                    end

                    obj.chooseStringAndUpdate(i);
                    
                    if k == 1
                        %Removal Stage
                        removedEdgeIndices(s) = i;
                    else
                        % Adding Stage
                        addedEdgeIndices(s) = i;
                    end
                end
            end
            
            if (numel(intersect(addedEdgeIndices, removedEdgeIndices)) == iterativeStepSize)
                fprintf('INFO: Breaking iterative step due to equal removed and added strings\n');
                break;
            end
            
            % Remove edges as long as there is an improvement
            obj.setRemovalMode(true);
            
            condition = true;
            minL2 = obj.diffToBlankSquaredErrorSum;
            
            fprintf('Try to improve by Removal...\n');
            while condition
                [val, i] = obj.findBestString();
                if val < minL2
                    obj.chooseStringAndUpdate(i);
                    minL2 = obj.diffToBlankSquaredErrorSum;
                else
                    condition = false;
                end
            end
            
            % Add edges as long as there is an improvement
            obj.setRemovalMode(false);
            
            condition = true;
            minL2 = obj.diffToBlankSquaredErrorSum;
            
            fprintf('Try to improve by Addition...\n');
            while condition
                [val, i] = obj.findBestString();
                
                if val < minL2
                    obj.chooseStringAndUpdate(i);
                    minL2 = obj.diffToBlankSquaredErrorSum;
                else
                    condition = false;
                end
            end
            
            rmse = obj.getRmseValue();
                    
            if (rmse < bestRmseValue)
                bestRmseValue = rmse;
                numBadRuns = 0;
                pickedEdgesSequence = obj.pickedEdgesSequence;
            else
                numBadRuns = numBadRuns + 1;
            end
        end        
        fprintf('INFO: Iterative greedy could improve RMSE...\nFrom %16.16f\nTo %16.16f\n', pureGreedyBestRmse, bestRmseValue);
    end

	x = 0 * obj.x;
	x(pickedEdgesSequence) = 1;
end

function removeOvershoot(object, numEdges)
    object.setRemovalMode(true);
    for k = 1 : numEdges
        fprintf('Removing string %i of %i\n', k, numEdges);
        i = object.pickedEdgesSequence(end);
        
        fprintf('\tF1 when picking edge Nr. %i: %16.16f\n', i, object.FRemoving(i));
        object.chooseStringAndUpdate(i);
    end
end