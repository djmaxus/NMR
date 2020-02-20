function [P] = greedy(basis)
%GREEDY Greedy algorithm (Chaturantabut et al. 2010)
%   Detailed explanation goes here
[~,i] = max(abs(basis(:,1)));
P = sparse(size(basis,1),size(basis,2)); 
P(i,1) = 1;
wb = waitbar(0,'Greedy sampling');
for j = 2:size(basis,2)
    V = basis(:,1:j-1);
    c = (P'*V)\(P'*basis(:,j));
    r = basis(:,j) - V*c;
    [~,i] = max(abs(r));
    dP = zeros(size(r)); 
    dP(i) = 1;
    P(:,j) = dP;
    progress =(j-1)/(size(basis,2)-1);
    waitbar(progress,wb);
end
waitbar(1,wb,'Greedy sampling: done');
% waitfor(wb);
delete(wb);
end

