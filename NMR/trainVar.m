Ns = size(vals,1);
data = cell(Ns,1);
parfor i = 1:Ns
    bc_i = bc;
    bc_i.val = vals(i,:);
    [data{i}] = solveDif(Nx,bc_i,D,1);   
    fprintf([num2str(i),' ']);
end
fprintf('\n');
data_x = [];
for i = 1:Ns
    data_x = [data_x,data{i}.X];
end
