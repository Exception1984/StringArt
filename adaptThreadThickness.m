function ttAdapted = adaptThreadThickness(tt, fd)

    theorRes = fd / tt;
    res = pow2(ceil(log2(theorRes)));
    
    ttAdapted = fd / res;
end