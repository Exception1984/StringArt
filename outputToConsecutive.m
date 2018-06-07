function outputToConsecutive(x, resultFilePath, matrixPath)

    if nargin == 0
        matrixPath = 'D:\Development\StringArt_Public_SourceCode\data\matrix_4096_256Pins';
        outputPath = 'D:\Development\StringArt_Public_SourceCode\output\einstein\einstein_NH-256_DW-512_WS-8.txt';
        resultFilePath = 'D:\Development\StringArt_Public_SourceCode\output\einstein\einstein_NH-256_DW-512_WS-8_RESULT.txt';
        
        % Read x
        fileID = fopen(outputPath, 'r');
        x = double(cell2mat(textscan(fileID, '%d')));
        fclose(fileID);
        
        temp = zeros(130560, 1);
        temp(x) = true;
        x = temp;
    end

    [~, ~, stringList, ~, ~, outsideMask] = findConsecutivePath(x, matrixPath);
    writeResultFile(stringList, outsideMask, resultFilePath);
end