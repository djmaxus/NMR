startup;
%%
Nx = 1000;
bc = genBC(1,[0,0]);
k = 1;
model = @(u) [1+u^2*k,2*u*k];
D  = genDiffusion(Nx,0.9,1,0.5e-1,model);
q = genSource(Nx,0.5,0);
generate = @(p) genProblem('TRAIN',Nx,genParam(D,bc,q),genInit(Nx,p(1),p(2),p(3)));
%%
T = 1;
Nt = 100;
dt = T/Nt;
newton_param = struct(...
    'iters'    , 10,...
    'iters_min', 3,...
    'tol'      , 1e-3,...
    'tol_rel'  , 1e-3,...
    'tol_g'    , 1e-3,...
    'tol_rel_g', 1e-3 ...
    );
nlin_newton = @(x,Res,Jac) newton(x,Res,Jac,newton_param);
%% Gen dataset
i0 = lhsdesign(10,3);
for i = 1:size(x0,2)
    problem_ = generate(x0(:,i));
    problem(i) = solveDif(T,Nt,problem_,nlin_newton);
end
