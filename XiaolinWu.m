%Matlab optmized version of XiolinWu line algorithm. No loops.
%Format:
%               [x y]=XiolinWu(x1,y1,x2,y2)
%
%Input:
%               (x1,y1): Start position
%               (x2,y2): End position
%
%Output:
%               x y: the line coordinates from (x1,y1) to (x2,y2)
%               c  : intensity of the line
%
%Usage example:
%               [x y]=bham(1,1, 10,-5);
%               plot(x,y,'or');

function [x, y, c] = XiaolinWu(x1, y1, x2, y2)
    dx = x2 - x1;
    dy = y2 - y1;
    %preallocate memory for x,y, and c
    x=zeros(floor(2*sqrt(dx^2+dy^2)),1);
    y=zeros(size(x));
    c=zeros(size(x));
    
    swapped=false;
    if abs(dx) < abs(dy)         
      [x1, y1]=swap (x1, y1);
      [x2, y2]=swap (x2, y2);
      [dx, dy]=swap (dx, dy);
      swapped=true;
    end 
    if x2 < x1
      [x1, x2]=swap (x1, x2);
      [y1, y2]=swap (y1, y2);
    end 
    gradient = dy / dx;
    
    % handle first endpoint
    xend = round(x1);
    yend = y1 + gradient * (xend - x1);
    xgap = rfpart(x1 + 0.5);
    xpxl1 = xend;  % this will be used in the main loop
    ypxl1 = ipart(yend);
    x(1)=xpxl1; y(1)=ypxl1; c(1)=rfpart(yend) * xgap;
    x(2)=xpxl1; y(2)=ypxl1 + 1; c(2)=fpart(yend) * xgap;
    intery = yend + gradient; % first y-intersection for the main loop
    
    % handle second endpoint
    xend = round (x2);
    yend = y2 + gradient * (xend - x2);
    xgap = fpart(x2 + 0.5);
    xpxl2 = xend;  % this will be used in the main loop
    ypxl2 = ipart (yend);
    x(3)=xpxl2; y(3)=ypxl2; c(3)=rfpart (yend) * xgap;
    x(4)=xpxl2; y(4)=ypxl2 + 1; c(4)=fpart (yend) * xgap;
    
    % main loop
    k=5;
    for i =(xpxl1 + 1):(xpxl2 - 1)
        x(k)=i; y(k)=ipart (intery); c(k)=rfpart (intery);
        k=k+1;
        x(k)=i; y(k)=ipart (intery) + 1; c(k)= fpart (intery);
        intery = intery + gradient; 
        k=k+1;
    end
    
    %truncate the vectors to proper sizes
    x=x(1:k-1);
    y=y(1:k-1);
    c=c(1:k-1);    
    if swapped         
      [x, y]=swap (x, y);
    end 

%integer part
function i=ipart(x)
    if x>0
        i=floor(x);
    else
        i=ceil(x);
    end

function r=round(x) 
    r= ipart(x + 0.5);

%fractional part
function f=fpart(x) 
    f=x-ipart(x);

function rf=rfpart(x) 
    rf= 1 - fpart(x);
    
function [x, y]=swap(x,y)
    a=x;
    x=y;
y=a;