function [alpha,delta]=sastailfit(x,L)

lg_x=log(abs(x));
nbins=50;
if ~isempty(L)
    L=log(L);
    edges=L(1):(L(2)-L(1))/nbins:L(2);
else
    qnt=(quantile(lg_x,[0.8, 0.95]));
    edges=qnt(1):(qnt(2)-qnt(1))/nbins:qnt(2);
end

bins=(edges(1:end-1)+edges(2:end))/2;
lg_epdf=histcounts(lg_x,edges,'normalization','pdf');

X=[bins.', ones(nbins,1)];

coeff=(X.'*X)\(X.'*lg_epdf);

grad= coeff(1);
c= coeff(2);

alpha=-grad-1;
Ca=gamma(alpha)*sin(alpha*pi/2)/pi;
delta=(exp(c)/(2*alpha*Ca)).^(1/alpha);
end