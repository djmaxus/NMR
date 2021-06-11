%% set problem
% dimensionality
Nx = 1000;
% BC
bc = genBC(1,[0,0]);
% diffusion
D = genDiffusion(Nx,0.9,1,0.5e-1,@(u)([1,0]));
% source term
q = genSource(Nx,0.6,1e+3);
% all params
param = genParam(D,bc,q);
% initial state
x0 = genInit(Nx,0.5,2,1);
% whole problem
problem = genProblem('FOM',Nx,param,x0);

%% set nonlinear solver
newton_param = struct(...
    'iters'    , 10,...
    'iters_min', 3,...
    'tol'      , 1e-3,...
    'tol_rel'  , 1e-3,...
    'tol_g'    , 1e-3,...
    'tol_rel_g', 1e-3 ...
    );

nlin_solver = @(x,Res,Jac) newton(x,Res,Jac,newton_param);

%% ODE settings
T = 1;
Nt = 100;
dt = T/Nt;

%% Solve problem numerically
[problem.sol] = solveDif(T,Nt,[dt,dt],problem,nlin_solver);
close all;
plotSol(problem.sol, figure);

%% Analytical solution
% generate linear matrix
J = buildJacobian(zeros(Nx,2),1,problem.param);
% u_t + Jx * u = q
Jx = full(J.Jx); % symmetric
% spectrum
[U,S,~] = svd(Jx); S = diag(S);
% inverse
S_inv = 1./S;
L_inv = U * diag(S_inv) * U';
Lp_inv = diag(1./S) * U';
% exact solution
thres = 1e-6; % sparse threshold
pos = problem.param.q.pos;
qval = problem.param.q.val(pos);
t = problem.sol.t;
a0 = U' * x0;
a0(abs(a0) < thres) = 0;
exact = problem.sol;
% forced term
p = U';
p = p(:,pos);
p(abs(p) < thres) = 0;
p = -qval * p;
dX = zeros(size(exact.X));
for n = 1:length(t)
    time = exp(-S*t(n));
    a = a0 .* time + p .* S_inv .* (1-time);
    dX(:,n) = U * a;
end
exact.X = dX;
plotSol(exact,figure);
plotComp(exact, problem.sol, D, figure());
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
