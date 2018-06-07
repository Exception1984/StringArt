function b = buildImageVector(img, domainWidthPixels)
    assert(size(img, 1) == size(img, 2));
    assert(size(img, 1) <= domainWidthPixels);

    % Take the transpose of the image to retrieve the pixels row by row
    % when applying the (:) operator
    imgTrans = img';

    imVector = imgTrans(:);
    
    v = imVector;
    m = domainWidthPixels * domainWidthPixels;
    n = 1;

    row_ind = 1 : m;
    col_ind = ones(size(row_ind, 1), 1);
    
    b = sparse(row_ind, col_ind, v, m, n);
end