%% figures
fig_def_prop = {'WindowStyle','docked'};
new_figure = @(name)(figure('Name',name,fig_def_prop{1:end}));
%% solver
bc.isKindI = @(kind)(-kind + 2);
q_init = @(pos,val)(struct('pos',pos,'val',val));
D_init = @(amp,period,scale)(struct('A',amp,'P',period,'S',scale));
