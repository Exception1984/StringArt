function [x, y, v] = drawLine(x1, y1, x2, y2, xOff, yOff, numLine, numLineTotal, domainMin, domainMax)
    xOffset = 0;
    yOffset = 0;
    
    if nargin > 4
        xOffset = xOff;
        if nargin > 5
            yOffset = yOff;
            if nargin > 6
                fprintf('Line Nr. %d of %d\n', numLine, numLineTotal);
            end
        end
    end

    %[x, y, v] = antialiasedLine(x1, y1, x2, y2);
    %[x, y, v] = bresenham(x1, y1, x2, y2);
    [x, y, v] = XiaolinWu(x1 + xOffset, y1 + yOffset, x2 + xOffset, y2 + yOffset);
    
    x = x - xOffset;
    y = y - yOffset;

    mask = (x >= domainMin) & (x <= domainMax) & (y >= domainMin) & (y <= domainMax);
    mask = mask & (v > 0);
    
    % Test cast to int32
    x = int32(x(mask));
    y = int32(y(mask));
    v = v(mask);
    
    if size(x, 2) == 0
        x = x';
        y = y';
        v = v';
    end
end