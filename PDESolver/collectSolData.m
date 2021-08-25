function [X,W,Jw,Y] = collectSolData(sol)
    X = sol.X;

    W = [sol.nlin.residual];
    W = [W.W];

    Jw = [sol.nlin.jacobian];
    Jw = [Jw.Jf];
    
    nlin = sol.nlin(1);
    Y = nlin.residual(1).W;
    for nlin = sol.nlin
        Y = [Y,nlin.residual(end).W];
    end
end
