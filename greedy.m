function [idx, P] = greedy(basis,n)
    idx = [];
    P = [];
    n = max(n,size(basis,2));
    N = size(basis,1);
    Nr = size(basis,2);
    for i = 1:Nr
        if i < mod(n,Nr) + 1
            ni = floor(n/Nr) + 1;
        else
            ni = floor(n/Nr);
        end
        if i == 1
            er = basis(:,1);
        else
            M = P' * basis(:,1:i-1);
            [u,s,v] = svd(M,'econ');
            M = v * inv(s) * u';
            M = eye(N) - basis(:,1:i-1) * M * P';
            er = M * basis(:,i);
        end
        for j = 1:ni
            cost = abs(er);
            [~,idx_cur] = sort(cost,'descend');
            for k = idx_cur'
                k_found = false;
                for p = idx
                    if k == p
                        k_found = true;
                        break;
                    end
                end
                if ~k_found
                    break;
                end
            end
            idx = [idx,k];
            dP = zeros(N,1);
            dP(k) = 1;
            P = [P,dP];
        end
    end
end

function [P] = greedy_old(basis)
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
