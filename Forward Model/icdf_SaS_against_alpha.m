clc; clear; close all;

% Finds the 'perc' quantile of SaS distributions, with varying alpha 

alpha=1.1:0.01:2;
icdfsas=zeros(1,length(alpha));
perc=0.99;

for i=1:length(alpha)
    pd=makedist('Stable','alpha',alpha(i),'beta',0,'gam',2,'delta',0);
    icdfsas(i)=icdf(pd,perc);
end
figure
plot(alpha,icdfsas);
xlabel('alpha')
ylabel(['icdf(',num2str(perc),')'])