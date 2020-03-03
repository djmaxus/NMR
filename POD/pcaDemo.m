clear all; close all; clc;
xC = [2;1;];
sig = [2;0.5];
theta = pi/3;
R = [ cos(theta),-sin(theta); 
      sin(theta), cos(theta) ];
  
nPoints = 2000;
X = R * diag(sig) * randn(2,nPoints) + diag(xC) * ones(2,nPoints);
xmin = min(X,[],'all');
xmax = max(X,[],'all');

f = figure;
f.WindowState = 'maximized';

subplot(1,2,1);
scatter(X(1,:),X(2,:),'k.','LineWidth',2);
hold on;
box on;
grid on;
axis square;
axis([xmin xmax xmin xmax]);

waitforbuttonpress;

%%
Xavg = mean(X,2);
B = X - Xavg;

[U,S,V] = svd(B/sqrt(nPoints-1),'econ');

theta = (0:0.01:1.0)*2*pi;
L = U*S;
[Xstd] = U*S*[cos(theta);sin(theta)];

for j = 1:3
plot(Xavg(1) + j*Xstd(1,:),Xavg(2) + j*Xstd(2,:),'r-','LineWidth',1.5);
end
plot([Xavg(1),Xavg(1) + L(1,1) ],[Xavg(2),Xavg(2) + L(2,1)],'-+g','LineWidth',2.5);
plot([Xavg(1),Xavg(1) + L(1,2) ],[Xavg(2),Xavg(2) + L(2,2)],'-+g','LineWidth',2.5);

waitforbuttonpress;
%% outliers
Y = R' * diag([0.1;0.1]) * randn(2,100) + diag([-4;6]) * ones(2,100);
X = [X,Y];
subplot(1,2,2);
scatter(X(1,:),X(2,:),'k.','LineWidth',2);
hold on;
box on;
grid on;
axis square;
xmin = min(X,[],'all');
xmax = max(X,[],'all');
axis([xmin xmax xmin xmax]);

waitforbuttonpress;
%% 
Xavg = mean(X,2);
B = X - Xavg;

[U,S,V] = svd(B/sqrt(size(X,2)-1),'econ');

L = U*S;
[Xstd] = U*S*[cos(theta);sin(theta)];

for j = 1:3
plot(Xavg(1) + j*Xstd(1,:),Xavg(2) + j*Xstd(2,:),'r-','LineWidth',1.5);
end
plot([Xavg(1),Xavg(1) + L(1,1) ],[Xavg(2),Xavg(2) + L(2,1)],'-+g','LineWidth',2.5);
plot([Xavg(1),Xavg(1) + L(1,2) ],[Xavg(2),Xavg(2) + L(2,2)],'-+g','LineWidth',2.5);

drawnow;

%%
% 
% f = figure;
% f.WindowState = 'maximized';
% scatter(X(1,:),X(2,:),'k.','LineWidth',2);
% hold on;
% box on;
% grid on;
% axis square;
% xmin = min(X,[],'all');
% xmax = max(X,[],'all');
% axis([xmin xmax xmin xmax]);
% drawnow;
% %% 
% 
% Xavg = mean(X,2);
% B = X - Xavg;
% 
% [Q,B] = l1pca(B/sqrt(size(X,2)-1),2);
% 
% theta = (0:0.01:1.0)*2*pi;
% L = Q*S;
% 
% plot([Xavg(1),Xavg(1) + L(1,1) ],[Xavg(2),Xavg(2) + L(2,1)*S(1,1)],'-+g','LineWidth',1.5);
% plot([Xavg(1),Xavg(1) + L(1,2) ],[Xavg(2),Xavg(2) + L(2,2)*S(2,2)],'-+g','LineWidth',1.5);
% 
% drawnow;
