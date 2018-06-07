function matrix = matrixFromFile(filePath, dataType)

    type = 'double';

    if nargin > 1
        type = dataType;
    end
    
    % Read matrix parameters
    fileID = fopen([filePath '.txt']);
    params = textscan(fileID, '%s %s %d');
    params = params{3};
    m = double(params(1));
    n = double(params(2));
    fclose(fileID);

    % Load matrix
    fileID = fopen([filePath '.bin']);
    matrix = fread(fileID,[m n], type);
    fclose(fileID);
    
    %Cast Matrix
    switch type
        case 'double'
            matrix = double(matrix);
        case 'int32'
            matrix = int32(matrix);
        case 'uint32'
            matrix = uint32(matrix);
            
    end
end