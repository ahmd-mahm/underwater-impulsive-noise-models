function [h,dB_bins]=cdfquant(x,nbins,L)

% Computes and plots the empirical CDF of the data vector x.
%
% Inputs: 
% x : input vector of observations nbins: number of bins in the CDF
% nbins : number of bins in the computed CDF
% L : (optional argument) A 2-dimensionl vector, where L(1)and L(2) are two such points on the
% CDF of x that determine the elements of x from which the empirical CDF is
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
    qnt=quantile(x,L);
else
    qnt=quantile(x,[0, 0.99]);
end
edges=qnt(1):(qnt(2)-qnt(1))/nbins:qnt(2);
dB_bins=(edges(1:end-1)+edges(2:end))/2;

h=histogram(x,edges,'normalization','cdf');
grid on

% figure
% histogram(x,edges,'normalization','pdf')
% set(gca,'xscale','log');
% set(gca,'yscale','log');
end