function [W] = buildFlowRes(X,D,bc,varargin)
    N = size(X,1);
    
    if length(varargin) == 1
        nfa = varargin{1};
        mask = nfa.mask;
        approx = @(f)(nfa.eval(f));
    else
        mask = 1:N+1;
        approx = @(f)(f);
    end
    
    
    % calculation
    W = zeros(N+1,1);
    for i = mask
        W(i) = calcFlowRes(X(:,2),D,bc,i);
    end    
    % gappy
    W = approx(W);
end
