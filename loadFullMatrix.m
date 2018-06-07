function A = loadFullMatrix(matrixPath, imSize)

    
    % Read parameters of filtered A
    fileID = fopen([matrixPath '/A_s=' num2str(imSize) '.txt']);
    params = textscan(fileID, '%s %s %d');
    params = params{3};
    m = double(params(1));
    n = double(params(2));
    nnz = double(params(3));
    fclose(fileID);

    % Load filtered version of A
    fileName = [matrixPath '/A_s=' num2str(imSize) '.bin'];
    if exist(fileName, 'file') == 2
        fileID = fopen([matrixPath '/A_s=' num2str(imSize) '.bin']);
        data = fread(fileID,[nnz 3],'double');
        fclose(fileID);

        A = sparse(data(:,1), data(:,2), data(:,3), m, n);
    else
        fileID = fopen([matrixPath '/A_s=' num2str(imSize) '_i.bin']);
        i = fread(fileID, [nnz 1], 'uint32');
        fclose(fileID);

        fileID = fopen([matrixPath '/A_s=' num2str(imSize) '_j.bin']);
        j = fread(fileID,[nnz 1], 'uint32');
        fclose(fileID);

        fileID = fopen([matrixPath '/A_s=' num2str(imSize) '_v.bin']);
        v = fread(fileID,[nnz 1], 'double');
        fclose(fileID);

        A = sparse(i, j, v, m, n);
    end
end