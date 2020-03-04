% II: residual and Jacobian snapshots and POD
data_r = [];
data_j = [];
for i = 1:size(vals,1)
    bc.val = vals(i,:);
    [data] = solveDif(Nx,bc,D,2,POD);
    data_r = [data_r,data.Rx];
    data_j = [data_j,data.Jnz];
end
