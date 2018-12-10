clc; clear; close all

% Populates uniformly distributed points in a circle
% Plots the output figure too

rho=10;     % radius of circle
N=10000;    % number of points
rho_p=sqrt(rand(1,N))*rho;
theta_p=rand(1,N)*2*pi;

figure
theta=0:2*pi/500:2*pi;
[x,y]=pol2cart(theta,repmat(rho,1,length(theta)));
plot(x,y,'-k','linewidth',2)
hold on

[xp,yp]=pol2cart(theta_p,rho_p);
plot(xp,yp,'.','markersize',2)
grid on
axis equal