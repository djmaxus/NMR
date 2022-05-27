function [im] = plotField(ax, Field, T, X, name, title_settings)
%PLOTFIELD Summary of this function goes here
%   Detailed explanation goes here

cla(ax);

im = imagesc(ax,'CData',Field,'XData',T,'YData', X);

colorbar(ax);
colormap(ax,'jet');

axis(ax,'equal');

xlim(ax,[0, 1]);
ylim(ax,[0, 1]);

label_settings = {'interpreter','latex'};
xlabel(ax,'$t$',label_settings{:});
ylabel(ax,'$x$',label_settings{:});

title(ax, name, title_settings{:});
end

