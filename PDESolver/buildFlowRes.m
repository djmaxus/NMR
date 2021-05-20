function [W] = buildFlowRes(X,D,bc,POD)
    N = size(X,1);
    nmr = ~isempty(POD);

    % mask for gappy-POD
    mask = 1:N+1;
    if nmr
        mask = POD.idx_mask;
    end
    
    % calculation
    W = zeros(N+1,1);
    for i = mask
        W(i) = calcFlowRes(X(:,2),D,bc,i);
    end
    
    % gappy
    if nmr
        W = applyOblique(W,POD);
    end 
end
