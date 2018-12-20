clc; clear all; close all

% Ascertains the derived log-normal distribution of I_t with a simulated distribution of the latter. 

N=10^7;
I_t_mean=40;
I_t_var=I_t_mean*0.1;
I_t=randn(N,1)*sqrt(I_t_var)+I_t_mean; % in dB

I_t=10.^(I_t/10);

figure
nbins=100;
histogram(I_t,nbins,'normalization','pdf')
xlabel('i_T')
ylabel('pdf')
grid on
hold on

x=0:0.1:max(I_t);
I_t_theory=pdf('lognormal',x,I_t_mean/(10*log10(exp(1))),sqrt(I_t_var)/(10*log10(exp(1))));
plot(x,I_t_theory,'-k','linewidth',2)