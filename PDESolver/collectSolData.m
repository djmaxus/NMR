function [X,W,Jw] = collectSolData(sol)
    X = sol.X;

    W = [sol.nlin.residual];
    W = [W.W];

    Jw = [sol.nlin.jacobian];
    Jw = [Jw.Jf];
end
