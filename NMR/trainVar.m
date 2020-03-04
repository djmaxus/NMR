data_x = [];
for i = 1:size(vals,1)
    bc.val = vals(i,:);
    [data] = solveDif(Nx,bc,D,1);
    data_x = [data_x,data.X];
end



