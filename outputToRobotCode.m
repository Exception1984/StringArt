clc;
clear;

matrixPath = 'D:\Development\StringArt_Public_SourceCode\data\matrix_4096_256Pins';
outputPath = 'D:\Development\StringArt_Public_SourceCode\output\einstein\einstein_NH-256_DW-512_WS-8.txt';
    
% Read x
fileID = fopen(outputPath, 'r');
x = double(cell2mat(textscan(fileID, '%d')));
fclose(fileID);

temp = zeros(130560, 1);
temp(x) = true;
x = temp;

[x, edgeCodes, stringList, hookSequence, hookApproachLeftSide, outsideMask] = findConsecutivePath(x, matrixPath);

%Write Result File
fileID = fopen(filePath, 'w');

% Write number of hooks
fprintf(fileID, '256\t0\t0\n');

% Write Hook Positions
hookPos = getUsedHookPositions();
for i = 1 : size(hookPos, 1)
    fprintf(fileID, '%d\t%16.16f\t%16.16f\n', i, hookPos(i, 1), hookPos(i, 2));
end

% Write x vector
nze = find(x);
for i = 1 : numel(nze)
    fprintf(fileID, '%d\t%d\t%d\n', edgeCodes(nze(i), 1), edgeCodes(nze(i), 2), x(nze(i)));
end

% Write number of edges
fprintf(fileID, '%d\t0\t0\n', sum(x));

dirs = [0; 1];
sides = [1; 0];

%hookLeft = double(hookApproachLeftSide) + 1;
outMask = double(outsideMask) + 1;
for i = 1 : size(stringList, 1)
    fprintf(fileID, '%d\t%d\t%d\t%d\t%d\t%d\n', 0, stringList(i, 1), stringList(i, 2), stringList(i, 3), stringList(i, 4), sides(outsideMask(i)));
end

fclose(fileID);