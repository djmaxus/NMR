function [out] = calcMSE(U,S,data)
out = zeros(1,size(S,2));
message = 'MSE calculation';
wb = waitbar(0,message);
for n = 1:size(S,2)-1
    waitbar(n/(size(S,2)-1),wb);
    A = U(:,1:n);
    restore = A * A' * data;
    delta = data - restore;
    delta = sum(delta.^2,'all')/numel(data);
    out(n) = delta;
end
close(wb);
end