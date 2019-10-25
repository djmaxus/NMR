addpath('./ResSim');
startMRSTmodules;
Nx = 21;
perm = 3.2 * ones(Nx^2,1);
obj = FOM(perm);
obj = obj.makeSnapshots();
%% prepare dataset
data = obj.snapshots;
init = data(:,1);
data = data - init;
mu = mean(data,2);
data = data - mu;
Ncells = obj.G.cells.num;

Pmax = max(data(1:Ncells,:),[],'all');
Pmin = min(data(1:Ncells,:),[],'all');
Smax = max(data(Ncells+1:end,:),[],'all');
Smin = min(data(Ncells+1:end,:),[],'all');
data(1:Ncells,:) = (data(1:Ncells,:) - Pmin) ./ (Pmax - Pmin);
data(Ncells+1:end,:) = (data(Ncells+1:end,:) - Smin) ./ (Smax - Smin);
POD = makePOD(data,0.9);
out_max = calcMSE(POD.U,POD.S,data);
% out_min = out_max;
% for i = 1:size(data,2)
%     temp = calcMSE(POD.U,POD.S,data(:,i));
%     out_max = max(temp,out_max);
%     out_min = min(temp,out_min);
% end
% hold on;
% semilogy(out_min);
semilogy(out_max);

[U,S,V] = svd(data);
Sx = zeros(size(S'));
diagSx = diag(S);
diagSx = eye(length(diagSx))./ diagSx;
Sx(1:end,1:length(diagSx))= diagSx;
Dx = V * Sx * U';
