% II: residual and Jacobian snapshots and POD
Ns = size(vals,1);
data = cell(Ns,1);
parfor i = 1:Ns
    bc_i = bc;
    bc_i.val = vals(i,:);
    [data{i}] = solveDif(Nx,bc_i,D,2,POD);
    disp(i);
end
data_r = [];
data_j = [];
for i = 1:Ns
    data_r = [data_r,data{i}.Rx];
    data_j = [data_j,data{i}.Jnz];
end
