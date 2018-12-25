clc; clear; close all

% accept-reject approach for point picking on the unit circle, such that
% the received pressure samples IID SaS samples.

d=5;        % sensor depth
h=20;       % height of water column
N=10^4;     % number of points
alpha=10;   % absorption coefficient in dB/km
x_max=5000;

It_mean_dB=100; % in dB
It_var_dB=(It_mean_dB*0.1/3)^2;
It_dB=It_mean_dB+(randn(1,N))*sqrt(It_var_dB); % log-normal distribution of intensity

alp=1.5;    % alpha of the SaS dist.
perc=1-10^-12;  % percentile match for the SaS dist. corresponding to Ir_dB = It_mean_dB - 20log(r) - alpha*(r/1000)

 
pd=makedist('Stable','alpha',alp,'beta',0,'gam',1,'delta',0);
icdfsas=icdf(pd,perc);
Pt_mean=10.^(It_mean_dB/20);
Pr_mean = 10.^((It_mean_dB- 20*log(h-d) - alpha*((h-d)/1000))/20); % in dB: Ir_dB = It_dB - 20log(r) - alpha*(r/1000)
del=Pr_mean/icdfsas;
pd=makedist('Stable','alpha',alp,'beta',0,'gam',del,'delta',0);
disp(['alpha = ',num2str(alp)]);
disp(['delta = ',num2str(del)]);

Pr=random(pd,1,N);
Pr_abs=abs(Pr);
Ir_dB=20*log10(Pr_abs);


%syms Ir_dB It_dB r alpha
%eqn = Ir_dB == It_dB - 20*log(r) - alpha*(r/1000);
%r= solve(eqn,r)
if alpha~=0
    r=(20000*wrightOmega(It_dB/20 - Ir_dB/20 - log(20000/alpha)))/alpha;
else
    r=10.^(It_dB-Ir_dB)/20;
end
x=sqrt(r.^2-(h-d)^2);
phi=2*pi*rand(1,N);
x_phase=x.*exp(1i*phi);


% in dB : Ir_dB = It_dB - 20log(r) - alpha*(r/1000). => alpha is in dB/km
% linear: Ir = It * (r^-2) * 10^(- alpha*r/(1000*10))

% It=10.^(It_dB/10);
% Ir=It.*(r.^(-2)); % adds spreading loss
% Ir = Ir.*10.^(-alpha*r/(1000*10)); % adds absorption


ppickingcircle(x_phase,x_max)

nbins = 100;
[lg_epdf,~,lg_bins]=logloghistquant(Pr_abs,nbins,[10^-4,1-10^-4]);
figure
bar(lg_bins,lg_epdf,'basevalue',min(lg_epdf)*(1-0.2*sign(min(lg_epdf))));
%set(gca,'xscale','log');
%set(gca,'yscale','log');
hold on

[alp_est,del_est,qnt]=sastailfit(Pr_abs,[0.9,1-10^-3]);
disp(['alpha = ',num2str(alp_est)]);
disp(['delta = ',num2str(del_est)]);

pdd = makedist('Stable','alpha',alp_est,'beta',0,'gam',del_est,'delta',0);
est_pdf=pdf(pdd,exp(lg_bins));
pdf_act=pdf(pd,exp(lg_bins));
plot(lg_bins,log(2*est_pdf),lg_bins,log(2*pdf_act),'linewidth',2)
plot(log(qnt),log(2*pdf(pdd,qnt)),'xk','markersize',12,'linewidth',2)

%histogram(Pr,nbins,'normalization','pdf');


function ppickingcircle(x_cmp,rho)

figure
theta=0:2*pi/500:2*pi;
[x,y]=pol2cart(theta,repmat(rho,1,length(theta)));
plot(x,y,'-k','linewidth',2)
hold on

plot(real(x_cmp),imag(x_cmp),'.','markersize',4)
grid on
axis equal
end

function [lg_epdf,edges,lg_bins]=logloghistquant(x,nbins,L)
if nargin==3
    qnt=log(quantile(x,L));
else
    qnt=log(quantile(x,[0, 0.99]));
end
edges=qnt(1):(qnt(2)-qnt(1))/nbins:qnt(2);
lg_bins=(edges(1:end-1)+edges(2:end))/2;

edges=exp(edges);

lg_epdf=log(histcounts(x,edges,'normalization','pdf')).';


% figure
% histogram(x,edges,'normalization','pdf')
% set(gca,'xscale','log');
% set(gca,'yscale','log');
end