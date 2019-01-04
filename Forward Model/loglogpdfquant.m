function [h,bins]=loglogpdfquant(x,nbins,L)

% Computes and plots the empirical PDF of the observations in x on a log-log
% scale.
%
% Inputs: 
% x : input vector of observations nbins: number of bins in the PDF
% nbins : number of bins in the computed PDF
% L : (optional argument) A 2-dimensionl vector, where L(1)and L(2) are two such points on the
% CDF of x that determine the elements of x from which the empirical PDF is
% computed. Note that L(1) and L(2) lie within  within [0,1] and L(2)>L(1)
% select the datapoints in x. Default value L=[0, 0.99]
%
% Outputs: 
% h : histogram object associated with the empirical PDF of x
% bins: a nbins-dimensional vector of bin locations
%
% -----
% Ahmed Mahmood (2019)

if nargin==3
    qnt=log(quantile(x,L));
else
    qnt=log(quantile(x,[0, 0.99]));
end
edges=qnt(1):(qnt(2)-qnt(1))/nbins:qnt(2);
lg_bins=(edges(1:end-1)+edges(2:end))/2;
edges=exp(edges);
bins=exp(lg_bins);

h=histogram(x,edges,'normalization','pdf');
set(gca,'xscale','log');
set(gca,'yscale','log');
grid on
title('log-log axis')
% figure
% histogram(x,edges,'normalization','pdf')
% set(gca,'xscale','log');
% set(gca,'yscale','log');
end