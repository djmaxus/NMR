%% set problem
Nx = 1000;
% BC
bc = genBC(1,[-1,1]);
D = genDiffusion(Nx,0.9,1,0.5e-1);
q = genSource(Nx,0.6,1e+3);
param = genParam(D,bc,q);

x0 = genInit(Nx,0.5,2,1);

problem = genProblem('FOM',Nx,param,x0);

%% set nonlinear solver
newton_param = struct(...
    'kmax'    , 10,...
    'kmin'    , 3,...
    'crit_rel', 1e-3,...
    'crit_abs', 1e-3...
    );

nlin_solver = newton_param;
nlin_solver.eval = @(x,Res,Jac) newton_common(x,Res,Jac,newton_param);

%% ODE settings
T = 1;
Nt = 100;
dt = T/Nt;

%% Solve problem
[problem.sol] = solveDif(T,Nt,[dt,dt],problem,nlin_solver);
%%
close all;
plotSol(problem.sol, figure);
%% another problem
problem_1 = genProblem('FOM-1',Nx,genParam(D,bc,genSource(Nx,0.598,1e+3)),x0);
[problem_1.sol] = solveDif(T,Nt,[dt,dt],problem_1,nlin_solver);
%% 
close all;
plotComp(problem_1.sol,problem.sol,param.D,figure);

%% build aDMDc
[X1,X2] = mergeForDMD(problem.sol.X,problem_1.sol.X);
Qa = full(repmat(problem.param.q.val,1,length(problem.sol.t)));
Qb = full(repmat(problem_1.param.q.val,1,length(problem_1.sol.t)));
[Q1,Q2] = mergeForDMD(Qa,Qb);
[A,B,C] = aDMDc(X1,X2,Q1,Q2);
[Psi,Lambda] = eigs(A,size(X1,2));
Lambda = diag(Lambda);

B(:,599) = 0.5*(B(:,598)+B(:,600));
C(:,599) = 0.5*(C(:,598)+C(:,600));
%% approximate of problem_1

X = x0;
for i = 1:Nt
    x = A*X(:,end) + B*problem_1.param.q.val + C*problem_1.param.q.val;
    X = [X,x];
end
approx = problem_1.sol;
approx.X = X;
plotComp(approx,problem_1.sol,param.D,figure);
