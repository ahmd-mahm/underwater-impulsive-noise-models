clc; clear; close all

N=10^6;
L=[1-10^-2, 1-10^-4];
nbins=100;

pd = makedist('Stable','alpha',1.5,'beta',0,'gam',3,'delta',0);
x=random(pd,N,1);


qnt=log(quantile(abs(x),[0.05, 0.9999]));
edges=exp(qnt(1):(qnt(2)-qnt(1))/nbins:qnt(2));
histogram(abs(x),edges,'normalization','pdf');
set(gca,'yscale','log')
set(gca,'xscale','log')

[alp,del,qnt]=sastailfit(x,L);

pdd = makedist('Stable','alpha',alp,'beta',0,'gam',del,'delta',0);
hold on
bins=(edges(1:end-1)+edges(2:end))/2;

est_pdf=pdf(pdd,bins);
plot(bins,2*est_pdf,'linewidth',2)

plot(qnt,2*pdf(pdd,qnt),'xk','markersize',12,'linewidth',2)