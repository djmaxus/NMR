function out = maskMatrix(mask)
    out = spalloc(size(mask,1),sum(mask,1),sum(mask,1));
    pos = 1;
    for i = 1:size(mask,1)
        if mask(i) == 1
            out(i,pos) = 1;
            pos = pos + 1;
        end
    end
end