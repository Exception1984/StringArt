function matrixDirPath = precomputeMatrix(numPins, threadThickness, frameDiameter, pinSideLength, superSamplingWindowWidth, dataPath) 
    highRes = round(frameDiameter / threadThickness);
    lowRes = highRes / superSamplingWindowWidth;
    
    matrixDirPath = fullfile(dataPath, matrixDirName(highRes, numPins));

    if ~(exist(matrixDirPath, 'dir') == 7)
        mkdir(matrixDirPath);
    end

    rowIndexFilePath = fullfile(matrixDirPath, ['A_s=' num2str(highRes) '_i.bin']);
    colIndexFilePath = fullfile(matrixDirPath, ['A_s=' num2str(highRes) '_j.bin']);
    valFilePath =      fullfile(matrixDirPath, ['A_s=' num2str(highRes) '_v.bin']);
    
    % Compute and write matrix files if they do not exist
    if  ~(exist(rowIndexFilePath, 'file') == 2) || ...
        ~(exist(colIndexFilePath, 'file') == 2) || ...
        ~(exist(valFilePath, 'file') == 2)

        % Compute Matrix A~
        [pinPos, fabricatableMask, i, j, v, m, n] = buildArcAdjacencyMatrix(numPins, threadThickness, frameDiameter, pinSideLength);
        
        % Wrte Pins Positions
        matrixToFile(pinPos, fullfile(matrixDirPath, 'pinPositions'), 'double');
        
        % Write mask of fabricatable edges
        matrixToFile(fabricatableMask, fullfile(matrixDirPath, 'fabricatableEdgesMask'), 'logical');
        
        A = sparse(i, j, v, m, n);
        
        % High Resolution Matrix
        handleMatrixContent(matrixDirPath, highRes, A, i, j, v, m, n);
        clear i j v
        
        % Create filtered matrix
        [i, j, v, m, n] = computeFilteredMatrices(lowRes, A);
        clear A
        A = sparse(i, j, v, m, n);
        
        % Low Resolution Matrix
        handleMatrixContent(matrixDirPath, lowRes, A, i, j, v, m, n);
        
        clear A i j v m n        
    end
end

function handleMatrixContent(matrixDirPath, res, A, i, j, v, m, n)
	numberOfNonZeros = numel(i);
        
    rowIndexFilePath = fullfile(matrixDirPath, ['A_s=' num2str(res) '_i.bin']);
    colIndexFilePath = fullfile(matrixDirPath, ['A_s=' num2str(res) '_j.bin']);
    valFilePath =      fullfile(matrixDirPath, ['A_s=' num2str(res) '_v.bin']);
    infoFilePath =     fullfile(matrixDirPath, ['A_s=' num2str(res) '.txt']);
    
    % Write Matrix
    writeBinaryFile(rowIndexFilePath, i, 'uint32');
    writeBinaryFile(colIndexFilePath, j, 'uint32');
    writeBinaryFile(valFilePath, v, 'double');

    % Write Info File
    infoText =  ['m = ' num2str(m) '\n' ...
                 'n = ' num2str(n) '\n' ...
                 'nnz = ' num2str(numberOfNonZeros)];
    writeTextFile(infoFilePath, infoText);

    % Compute Index Files
    AEdgeIndicesToPixelCodes = cell(n, 2);
    for k = 1 : n
        [row, ~, val] = find(A(:, k));
        AEdgeIndicesToPixelCodes{k, 1} = uint32(row);
        AEdgeIndicesToPixelCodes{k, 2} = val;
    end

    edgePixelIndices = cell2mat(AEdgeIndicesToPixelCodes(:, 1));
    matrixToFile(edgePixelIndices, fullfile(matrixDirPath, ['A_s=' num2str(res) '_EPI']), 'uint32');
    matrixToFile(cell2mat(AEdgeIndicesToPixelCodes(:, 2)), fullfile(matrixDirPath, ['A_s=' num2str(res) '_EPV']));
    matrixToFile(uint32(repelem((1 : n)', cellfun('length', AEdgeIndicesToPixelCodes(:, 1)))), fullfile(matrixDirPath, ['A_s=' num2str(res) '_CEI']), 'uint32');

    numberOfNonZeros = numel(edgePixelIndices);
    m = size(edgePixelIndices, 1);
	n = max(edgePixelIndices);
    
    filePath = fullfile(matrixDirPath, ['A_s=' num2str(res) '_ITI']);

    % Write Binary Index File
    writeBinaryFile([filePath '.bin'], uint32(edgePixelIndices), 'uint32');

    % Write Info File for ITI Index File
    infoText =  ['m = ' num2str(m) '\n' ...
                 'n = ' num2str(n) '\n' ...
                 'nnz = ' num2str(numberOfNonZeros)];
    writeTextFile([filePath '.txt'], infoText);
end

function writeBinaryFile(path, content, datatype)
    fileID = fopen(path, 'w');
    fwrite(fileID, content, datatype);
    fclose(fileID);
end

function writeTextFile(path, content)
    fileID = fopen(path, 'w');
    fprintf(fileID, content);
    fclose(fileID);
end