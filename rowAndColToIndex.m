function index = rowAndColToIndex(row, col, imageWidth)
    index = imageWidth * (row - 1) + col;
end