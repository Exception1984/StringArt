function [rowSpan, colSpan] = rowAndColSpan(row, col, numRows, numCols)

    if nargin < 4
        numCols = numRows;
    end

    n = numel(row);
    
    offRows = floor(numRows / 2);
    offCols = floor(numCols / 2);
    
    rowMin = row - offRows;
    colMin = col - offCols;

    stepRows = numRows - 1;
    stepCols = numCols - 1;
    
    numEntries = numRows * numCols;
    
    row_shift = repmat((0 : stepRows), numCols, 1);
    row_shift = repmat(row_shift(:)', n, 1);
    rowSpan = repmat(rowMin, 1, numEntries) + row_shift;
    
    col_shift = repmat((0 : stepCols), n, numRows);
    colSpan = repmat(colMin, 1, numEntries) + col_shift;
end