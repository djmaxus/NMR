%% load dataset
load 'data/gauss60_5+';
%% normalize dataset
mu = mean(dataset,2);
dataset = dataset - mu;
%% perform SVD
[U,S,V] = svd(dataset);
%% set gap parameters
% unmarred fraction
p = 0.1;
% unmarred number
n_unmarred = round(p * size(dataset,1));
num_eig_keep = 30;
%% set up masks
i_unmarred = randi([1,size(dataset,1)],n_unmarred,size(dataset,2));
mask = makeMask(dataset,n_unmarred);
data_masked = dataset .* mask;
%% marred average
semilogy(diag(S));
mu_masked = sum(data_masked,2)./ (sum(mask,2) + eq(sum(mask,2),0));
data_masked = data_masked .* eq(mask,1) + eq(mask,0) .* mu_masked;
max_iter = 20;
for num_iter = 1:max_iter
    [Um,Sm,Vm] = svd(data_masked);
    hold on;
    semilogy(diag(Sm));
    %% least squares problem
    A =  Um(:,1:num_eig_keep);
    a = A'*data_masked;
    %%
    data_masked = data_masked .* eq(mask,1) + eq(mask,0) .* (A*a);
end
%% check if entry is non-unique in array
function [flag] = is_repeat(entry_index,array)
sum_left = sum(eq(array(entry_index),array(1:entry_index-1)));
sum_right = sum(eq(array(entry_index),array(entry_index+1:end)));
flag = (sum_left + sum_right) > 0;
end

%% calculate mask
function mask = makeMask(dataset,n_unmarred)
mask = zeros(size(dataset));
i_unmarred = randi([1,size(dataset,1)],n_unmarred,size(dataset,2));
for i = 1:size(dataset,2)
    for j = 2:n_unmarred
        while is_repeat(j,i_unmarred(:,i))
            val = randi([1,size(dataset,1)],1,1);
            i_unmarred(j,i) = val;
        end
        mask(i_unmarred(j,i),i) = 1;
    end
end
end
function out = maskMatrix(mask)
    out = zeros(sum(mask,1),size(mask,1));
    pos = 1;
    for i = 1:size(mask,1)
        if mask(i) == 1
            out(pos,i) = 1;
            pos = pos + 1;
        end
    end
end