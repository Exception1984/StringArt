function [rows, cols, data, m, n, bFiltered] = computeFilteredMatrices(sideLength, A, b)
    if (nargin > 2)
        assert(size(A, 1) == size(b, 1), 'Size of A and b are incompatible!');
    end

    numPixels = size(A, 1);
    imageSize = sqrt(numPixels);

    filteredSize = sideLength; %round(scale * imageSize);
    
    if (nargin > 2)
        % Filter b
        bFiltered = full(reshape(b, [imageSize, imageSize]));
        bFiltered = imresize(bFiltered, [filteredSize filteredSize]);
        bFiltered = bFiltered(:);
    end
    
    % Filter A
    m = filteredSize * filteredSize;
    n = size(A, 2);
    
    %[X, Y] = meshgrid(1 : filteredSize, 0 : filteredSize - 1);
    
    %figure;
    %hold on;
    
    edgeImage = zeros(imageSize, imageSize);
    
    rows = zeros(0,1);
    cols = rows;
    data = rows;
    
    for i = 1 : n
        edgeImage = 0 * edgeImage;
        
        [j, k, v] = find(A(:, i));
        
        edgeImage(j) = v;
        edgeImage = edgeImage';
        filteredEdgeImage = imresize(edgeImage, [filteredSize filteredSize], 'box');

        [j, k, v] = find(filteredEdgeImage);
        
        %imshow(edgeImage);
        %imshow(filteredEdgeImage); 
        
        numElems = size(j, 1);
        
        rowNumbers = (j - 1) * filteredSize + k;
        colNumbers = repmat(i, numElems, 1);
        filteredData = v;
        
        %Update data
        rows = [rows' rowNumbers']';
        cols = [cols' colNumbers']';
        data = [data' filteredData']';
        
        fprintf('Iteration %d of %d\n', i, n);
    end
    
    %AFiltered = sparse(rows, cols, data, m, n);
end