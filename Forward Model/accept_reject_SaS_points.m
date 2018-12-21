clc; clear; close all

% accept-reject approach for point picking on the unit circle, such that
% the received pressure samples IID SaS samples.

d=5;        % sensor depth
h=20;       % height of water column
%c=1500;
N=10^5;     % number of points
alpha=10;   % absorption coefficient in dB/km
x_max=5000;

nbins=100;


%*** via the direct method
x=sqrt(rand(1,N))*x_max; % Uniform point picking in a circle
phi=2*pi*rand(1,N);
[x1,x2]=pol2cart(phi.',x.');
x_phase=x.*exp(1i*phi);

It_mean_dB=100; % in dB
It_var_dB=(It_mean_dB*0.1/3)^2;
It_dB=It_mean_dB+(randn(1,N))*sqrt(It_var_dB); % log-normal distribution of intensity
r=sqrt(x.^2+(d-h).^2);
% in dB : Ir_dB = It_dB - 20log(r) - alpha*(r/1000). => alpha is in dB/km
% linear: Ir = It * (r^-2) * 10^(- alpha*r/(1000*10))

%It=10.^(It_dB/10);
% Ir=It.*(r.^(-2)); % adds spreading loss
% Ir = Ir.*10.^(-alpha*r/(1000*10)); % adds absorption

Ir_dB = It_dB - 20*log(r) - alpha*(r/1000);
Ir=10.^(Ir_dB/10);
Pr=sqrt(Ir);

ppickingcircle(x_phase,x_max)

[epdf,edges,bins]=logxhistquant(Pr,nbins,[0,1-10^-4]);
figure
bar(bins,epdf);
set(gca,'xscale','log');
set(gca,'yscale','log');
hold on

[alp,del,qnt]=sastailfit(x,[0.8,1-10^-3]);
pdd = makedist('Stable','alpha',alp,'beta',0,'gam',del,'delta',0);

est_pdf=pdf(pdd,bins);
plot(bins,2*est_pdf,'linewidth',2)
plot(qnt,2*pdf(pdd,qnt),'xk','markersize',12,'linewidth',2)

%histogram(Pr,nbins,'normalization','pdf');





function ppickingcircle(x_cmp,rho)

figure
theta=0:2*pi/500:2*pi;
[x,y]=pol2cart(theta,repmat(rho,1,length(theta)));
plot(x,y,'-k','linewidth',2)
hold on

plot(real(x_cmp),imag(x_cmp),'.','markersize',2)
grid on
axis equal
end

function [epdf,edges,bins]=logxhistquant(x,nbins,L)
if nargin==3
    qnt=log(quantile(x,L));
else
    qnt=log(quantile(x,[0, 0.99]));
end
edges=qnt(1):(qnt(2)-qnt(1))/nbins:qnt(2);

edges=exp(edges);
bins=(edges(1:end-1)+edges(2:end))/2;

epdf=(histcounts(x,edges,'normalization','pdf')).';


% figure
% histogram(x,edges,'normalization','pdf')
% set(gca,'xscale','log');
% set(gca,'yscale','log');
end