function [Jf] = buildFlowJac(X,D,bc,varargin)
    N = size(X,1);
    
    if length(varargin) >= 1
        nfa = varargin{1};
        mask = nfa.mask;
        approx = @(f)(nfa.eval(f));
    else
        mask = 1:N+1;
        approx = @(f)(f);
    end
    
    
    
    Jf = zeros(N+1,2);
    for i = mask
        Jf(i,:) = calcFlowJac(X(:,2),D,bc,i,1:2);
    end
    
    % gappy;
    Jf = approx(Jf);
end

