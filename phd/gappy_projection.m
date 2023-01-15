function P = gappy_projection(dim_full,mask)
len_mask = length(mask);
P = sparse(mask,1:len_mask,1,dim_full,len_mask);
end

