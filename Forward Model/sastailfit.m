function [alpha,delta,qnt]=sastailfit(x,L)
% fits the tail of the empirical pdf computed from 'x', where L(1) and L(2)
% highlight the minimum and maximum percentile of the data to be considered
% in the tail. Note that L(1)<L(2) and both lie within [0,1]

x=abs(x);
nbins=50;
if nargin==2
    qnt=log(quantile(x,L));
else
    qnt=log(quantile(x,[0.99, 0.999]));
end

edges=qnt(1):(qnt(2)-qnt(1))/nbins:qnt(2);
bins=(edges(1:end-1)+edges(2:end))/2;
edges=exp(edges);

lg_epdf=log(histcounts(x,edges,'normalization','pdf')).';

bins=bins(~isinf(lg_epdf));
lg_epdf=lg_epdf(~isinf(lg_epdf));

X=[bins.', ones(length(bins),1)];
coeff=(X.'*X)\(X.'*lg_epdf);

grad= coeff(1);
c= coeff(2);

% figure
% bar(bins,lg_epdf,'basevalue',min(lg_epdf)*(1-0.2*sign(min(lg_epdf))))
% hold on
% plot(bins,grad*bins+c)
% plot(qnt,[-10,-10],'xk','linewidth',2)
% pause

alpha=-grad-1;
Ca=gamma(alpha)*sin(alpha*pi/2)/pi;
delta=(exp(c)/(2*alpha*Ca)).^(1/alpha);
qnt=exp(qnt);
end