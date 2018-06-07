function mseValue = mse(x, y)
    mseValue = sum((x - y).^2) / numel(x);
end