function [row, col] = indexToRowAndCol(index, imageWidth)
    row = floor((index - 1) / imageWidth) + 1;
    col = mod((index - 1), imageWidth) + 1;
end