function im = stretchHistogram(image)
    %Clamp Image if it is out of range
    image = max(0.0, min(1.0, image));

    mi = min(min(image));
    range = max(max(image)) - mi;
    im = (image - mi) / range;
end