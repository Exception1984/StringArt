function matrixToFile(matrix, filePath, dataType)

    type = 'double';
    
    if nargin > 2
        type = dataType;
    end

    % Write Binary File
    fileID = fopen([filePath '.bin'], 'w');
    fwrite(fileID, matrix, type);
    fclose(fileID);

    % Write Info File
    fileID = fopen([filePath '.txt'], 'w');
    fprintf(fileID, ['m = ' num2str(size(matrix, 1)) '\n']);
    fprintf(fileID, ['n = ' num2str(size(matrix, 2)) '\n']);
    fprintf(fileID, ['nnz = ' num2str(nnz(matrix)) '\n']);
    fprintf(fileID, ['type = ' type]);
    fclose(fileID);
end