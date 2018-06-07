function [x, y, v] = antialiasedLine(x1, y1, x2, y2, lw)
    lineWidth = 1;
    
    if (nargin > 4)
        lineWidth = lw;
    end

    width = abs(x2 - x1) + 1 + 2 * lineWidth;
    height = abs(y2 - y1) + 1 + 2 * lineWidth;
    
    origin = [min(x1, x2) min(y1, y2)];
    
    I = zeros(height, width);
    
    points = [x1, y1, x2, y2];
    offset = [1 + lineWidth, 1 + lineWidth];
    
    points(:,1:2) = points(:,1:2) - origin + offset;
    points(:,3:4) = points(:,3:4) - origin + offset;
    
    line = insertShape(I,'Line',points,'Color','white','LineWidth',lineWidth,'SmoothEdges',true);
    [row,col,v] = find(line(:,:,1));
    
    x = col - offset(1) + origin(1);
    y = row - offset(2) + origin(2);
end