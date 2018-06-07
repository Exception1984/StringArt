function highResIndex = lowResToHighResIndex(lowResIndex, lowResWidth, highResWidth)
    assert(mod(highResWidth, lowResWidth) == 0);
    
    filterWidth = highResWidth / lowResWidth;

    [rowLowRes, colLowRes] = indexToRowAndCol(lowResIndex, lowResWidth);

    offset = -0.5 * filterWidth + 1;
    
    [rowHighResSpan, colHighResSpan] = rowAndColSpan(filterWidth * rowLowRes + offset, filterWidth * colLowRes + offset, filterWidth, filterWidth);
    
    highResIndex = rowAndColToIndex(rowHighResSpan, colHighResSpan, highResWidth);
end