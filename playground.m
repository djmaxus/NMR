

%% set problem
% dimensionality
Nx = 1000;
% BC
bc = genBC(1,[-0.5,-0.5]);
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
clc;
[problem.sol] = solveDif(T,Nt,problem,nlin_solver);
plotSol(problem.sol, fig); fig;

%% Analytical solution
% generate linear matrix
J = buildJacobian(zeros(Nx,2),1,problem.param);
% u_t + Jx * u = q
Jx = full([J.Jx,J.Jxb;zeros(2,Nx+2)]);
% evaluate
exact = problem.sol;
X = calcGreenODE(Jx,[x0;bc.state(:,2)],problem.sol.t,-[problem.param.q.val;0;0]);
exact.X = X(1:end-2,:);
plotComp(exact, problem.sol, D, fig); fig;

%% try DMDc
X = exact.X;
dX = diff(X,1,2)./dt;
Q = repmat(-full(problem.param.q.val),1,Nt);
[A,U,S] = dmds(X(:,2:end),dX-Q,23);
dmded = problem.sol;
dmded.X = calcGreenODE(-A,x0,problem.sol.t,-problem.param.q.val,23);
dmded.X = real(dmded.X);
plotComp(dmded, exact, D, fig); fig;

%% another problem
problem_1 = genProblem('FOM-1',Nx,genParam(D,bc,genSource(Nx,0.5,1e+3)),x0);
[problem_1.sol] = solveDif(T,Nt,[dt,dt],problem_1,nlin_solver);
plotComp(problem_1.sol,problem.sol,param.D,fig);fig;

%% study
clc;
fig = new_figure('Result');
disp('Study start');

%% common parameters for all problems
Nx = 1000;                                      % dimensionality
bc = genBC(1,[0,0]);                            % boundary
D  = genDiffusion(Nx,0.9,1,0.5e-1,@(u)([1,0])); % diffusion
x0 = genInit(Nx,0.5,2,1);                       % initial state

%% Problem Constructors
makeFOM = @(q,name) genProblem(name,Nx,genParam(D,bc,q),x0);
makeROM_I = @(name,q,encoder) genProblem(name,Nx,genParam(D,bc,q),x0,encoder);
makeROM_II = @(name,q,encoder,nfa) genProblem(name,Nx,genParam(D,bc,q),x0,encoder,nfa.res,nfa.jac);

%% ODE settings
T = 1;
Nt = 100;
dt = T/Nt;

%% nonlinear solvers
newton_param = struct(...
    'iters'    , 10,...
    'iters_min', 3,...
    'tol'      , 1e-3,...
    'tol_rel'  , 1e-3,...
    'tol_g'    , 1e-3,...
    'tol_rel_g', 1e-3 ...
    );

nlin_newton = @(x,Res,Jac) newton      (x,Res,Jac,newton_param);
nlin_gauss  = @(x,Res,Jac) newton_gauss(x,Res,Jac,newton_param);

%% train set
train.q = genSource(Nx,0.6,1e+3);                         % specific property
    %% FOM
    train.FOM = makeFOM(train.q,'train.FOM');             % setup
    train.FOM.sol = solveDif(T,Nt,train.FOM,nlin_newton); % solve
    plotSol(fig,train.FOM.sol);                           % plot
    
    %% ROM-I
    % Encoder
    train.pod = PCA(train.FOM.sol.X,0.99,'init',x0,'sval_deg',1);
    train.encoder = train.pod.genEncoder();
    % Problem
    train.ROM.I = makeROM_I('train.ROM-I',train.q,train.encoder);
    train.ROM.nlin = @(x,Res,Jac) newton_gauss(x,Res,Jac,newton_param,train.encoder.decjac);
    % Solution
    train.ROM.I.sol = solveDif(T,Nt,train.ROM.I,train.ROM.nlin); % solve
    train.ROM.I.sol.X = train.encoder.decode(train.ROM.I.sol.X); % decode solution
    plotComp(fig,train.FOM.sol,train.ROM.I.sol);                 % plot
    
    %% ROM-II
    % Collect data
    [train.W,train.Jw] = collectNonlinData(train.FOM.sol);
    % Make POD
    N1 = floor(train.pod.npc() * 1.5);
    train.gappy.R = PCA(train.W,N1,'init',zeros(Nx+1,1));
    train.gappy.J = PCA(train.Jw,N1,'init',zeros(Nx+1,1));
    % Set masks
    N2 =floor(N1*1.5);
    idx = unique([train.gappy.R.calcGreedy(N2),...
                  train.gappy.J.calcGreedy(N2)]);
              
    train.gappy.R = train.gappy.R.buildGappy(idx);
    train.gappy.J = train.gappy.J.buildGappy(idx);
    
    train.nfa.res = train.gappy.R.genNFA();
    train.nfa.jac = train.gappy.J.genNFA();
    
    train.ROM.II = makeROM_II('ROM-II',train.q,train.encoder,train.nfa);
    train.ROM.II.sol = solveDif(T,Nt,train.ROM.II,train.ROM.nlin); % solve
    train.ROM.II.sol.X = train.encoder.decode(train.ROM.II.sol.X); % decode solution
    plotComp(fig,train.FOM.sol,train.ROM.II.sol);  
%% test set
test.q = genSource(Nx,0.55,1e+3);                       % specific property
    % FOM
    test.FOM = makeFOM(test.q,'test.FOM');          % setup
    test.FOM.sol = solveDif(T,Nt,test.FOM,nlin_newton); % solve
    plotComp(fig, train.FOM.sol, test.FOM.sol);                         % plot
    % ROM-I
    % ROM-II
    
%% Hypothesis
    % hypothesis #1 : build ROM-1 via DMD
    % hypothesis #2 : build ROM-1 via SPT

%% Comparison
plotComp(fig,train.FOM.sol,test.FOM.sol);

%% Collect ROM-II data
function [W,Jw] = collectNonlinData(sol)
    W = [];
    Jw = [];
    for nlin = sol.nlin
        for res = nlin.residual
            W = [W,res.W];
        end
        for jac = nlin.jacobian
            Jw = [Jw,jac.Jf];
        end
    end
end
