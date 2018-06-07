function [im, mask] = maskImage(inputImage, r)

    rows = size(inputImage, 1);
    cols = size(inputImage, 2);

    imageRangeY = (1 : cols) - 0.5 * cols - 0.5;
    imageRangeX = (1 : rows) - 0.5 * rows - 0.5;

    [x, y] = meshgrid(imageRangeX, imageRangeY);
    
    mask = x.^2 + y.^2 > r.^2;
    
    inputImage(mask) = 0;
    im = inputImage;
end