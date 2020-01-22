addpath('./ResSim');
startMRSTmodules;
Nx = 11;
perm = 3.2 * ones(Nx^2,1);
obj = FOM(perm);
obj = obj.makeSnapshots();

%% make POD and estimate quality
POD = makePOD(data,0.9);
out_max = calcMSE(POD.U,POD.S,data);
semilogy(out_max);
