function [alpha,delta]=sastailfit(x,L)
% fits the tail of the empirical pdf computed from 'x', where L(1) and L(2)
% highlight the percentile data to be considered in the tail. Note that
% L(1)<L(2) and both lie within [0,1]

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

lg_epdf=log(histcounts(abs(x),edges,'normalization','pdf')).';

% figure
% plot(bins,lg_epdf)

X=[bins.', ones(nbins,1)];
coeff=(X.'*X)\(X.'*lg_epdf);

grad= coeff(1);
c= coeff(2);

%hold on
%plot(bins,grad*bins+c)
hold on
plot(exp(bins),exp(grad*bins+c),'linewidth',2)


alpha=-grad-1;
Ca=gamma(alpha)*sin(alpha*pi/2)/pi;
delta=(exp(c)/(2*alpha*Ca)).^(1/alpha);
end