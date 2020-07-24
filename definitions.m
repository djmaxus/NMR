%% figures
fig_def_prop = {'WindowStyle','docked'};
new_figure = @(name)(figure('Name',name,fig_def_prop{1:end}));
%% solver
bc.isKindI = @(kind)(-kind + 2);
q_init = @(pos,val)(struct('pos',pos,'val',val));
D_init = @(amp,period,scale)(struct('A',amp,'P',period,'S',scale));
%% RBF
RBF.Gauss = @(x,x0,eps)(exp(-(eps*abs(x-x0)).^2));
RBF.Quad.Multi = @(x,x0,eps)(sqrt(1  + ( eps * abs(x-x0).^2 ) ) );
RBF.Quad.Inverse = @(x,x0,eps)((1  + ( eps * abs(x-x0).^2 ) ).^(-1) );
RBF.Quad.InverseMult = @(x,x0,eps)(sqrt(1  + ( eps * abs(x-x0).^2 ) ).^(-1) );
RBF.PolyHarm.Std = @(x,x0,k)(1 + abs(x-x0).^k);
