addpath('./ResSim');
startMRSTmodules;
Nx = 11;
perm = 3.2 * ones(Nx^2,1);
obj = FOM(perm);
obj = obj.makeSnapshots();

%% prepare dataset
data = obj.snapshots;
% substract initial
init = data(:,1);
data = data - init;
% substract mean
mu = mean(data,2);
data = data - mu;
disp(mean(data,2));
%divide by variation
sigma = sqrt(var(data,0,2));
data = data.\sigma;

%% make POD and estimate quality
POD = makePOD(data,int64(11));
out_max = calcMSE(POD.U,POD.S,data);
semilogy(out_max);
