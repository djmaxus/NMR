function [ax, leg] = plotProblem(ax,sol)
%PLOTPROBLEM Summary of this function goes here
%   Detailed explanation goes here
cla(ax);

hold(ax, 'on');

yyaxis(ax, 'left');
interpreter_latex = {'interpreter','latex'};
yll = ylabel(ax, '$$u$$',interpreter_latex{:});
yll.Rotation = 0;
%% boundary
bounds = [0,1];
plot(ax, bounds, sol.param.bc.val,'-.');

%% initial
Nx = sol.Nx;
dx = 1/Nx;
cell_coords = linspace(dx/2,1-dx/2,Nx);

x0 = sol.x0(:,1);
plot(ax, cell_coords,x0,'-b');

%% source
q_pos = sol.param.q.pos;
q_val = sol.param.q.val(q_pos);
q_label = ['$$q = ',num2str(q_val),'$$',newline(),'$$x_q = ',num2str(q_pos/Nx),'$$'];
xQ = xline(ax,cell_coords(q_pos),'--', q_label, interpreter_latex{:});
xQ.LabelOrientation = 'horizontal';

%% kappa
kappa = sol.param.D.val;
yyaxis(ax, 'right');
ylr = ylabel(ax, '$$\kappa$$',interpreter_latex{:});
ylr.Rotation = 0;
ylim(ax,[0,max(kappa)]);

edge_coords = linspace(0,1,Nx+1);

area(ax, edge_coords, kappa,'FaceColor','#999','LineStyle','none','FaceAlpha',0.5);

%%
hold(ax, 'off');

xlim(bounds);
axis(ax,'square');

xlabel('$$x$$',interpreter_latex{:});

leg = legend({'$$u_\mathrm{B}$$', '$$u(x,0)$$', '$$q$$', '$$\kappa$$'},interpreter_latex{:});
leg.Location = 'southoutside';
leg.NumColumns = 4;
leg.Box = 'on';

title(ax, 'Problem parameters');

end

