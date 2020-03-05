function [idx, P] = greedy(basis,n)
    idx = zeros(1,n);
    P = spalloc(size(basis,1),n,n);
    n = max(n,size(basis,2));
    N = size(basis,1);
    Nr = size(basis,2);
    nt = 0;
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
            nt = nt+1;
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
            idx(nt) = k;
            P(k,nt) = 1;
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
