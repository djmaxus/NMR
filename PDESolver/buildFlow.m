function [Jf] = buildFlow(X,D,bc,POD)
    N = size(X,1);
    nmr = ~isempty(POD);

    % gappy mask
    mask = 1:N+1;
    if nmr
        mask = POD.idx_mask;
    end
    
    Jf = zeros(N+1,2);
    for i = mask
        Jf(i,:) = calcFlowJac(X(:,2),D,bc,i,1:2);
    end
    
    % gappy;
    if nmr
        Jf = applyOblique(Jf.*[1,-1],POD).*[1,-1];
    end
end
