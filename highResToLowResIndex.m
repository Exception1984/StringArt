function lowResIndex = highResToLowResIndex(highResIndex, highResWidth, lowResWidth)
    assert(mod(highResWidth, lowResWidth) == 0);
    
    filterWidth = highResWidth / lowResWidth;

    [rowHighRes, colHighRes] = indexToRowAndCol(highResIndex, highResWidth);
    
    lowResIndex = rowAndColToIndex( ...
        floor((rowHighRes - 1) / filterWidth) + 1, ...
        floor((colHighRes - 1) / filterWidth) + 1, ...
        lowResWidth);
end