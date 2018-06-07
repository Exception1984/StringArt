function map = multiSampleCorrespondenceMap(domainWidth, windowSize)
    m = domainWidth * windowSize;
    numPixels = m * m;    
    
    k = domainWidth;
    numCorrespondenceValues = k * k;    
    
    row_ind = reshape(1 : numCorrespondenceValues, [k, k])';
    row_ind = imresize(row_ind, [m, m], 'box')';
    row_ind = row_ind(:);
    
    col_ind = (1 : numPixels)';

    v = (1 / (windowSize^2)) * ones(numPixels, 1);

    map = sparse(row_ind, col_ind, v, numCorrespondenceValues, numPixels);
end