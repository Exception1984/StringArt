function fabEdgesMask = loadFabricatableEdgesMask(matrixPath)

    fileID = fopen(fullfile(matrixPath, 'fabricatableEdgesMask.txt'), 'r');
    params = textscan(fileID, '%s %s %d');
    params = params{3};
    m = double(params(1));
    n = double(params(2));
    fclose(fileID);

    % Load matrix
    fileID = fopen(fullfile(matrixPath, 'fabricatableEdgesMask.bin'), 'r');
    fabEdgesMask = logical(fread(fileID,[m n], 'logical'));
    fclose(fileID);
end