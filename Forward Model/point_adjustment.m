clc; clear; close all

% corrects in the MS sense, the divergence between the Stable PDF and
% Simulated PDF by adjusting the Poisson points on the sea floor

d=5;
h=10;
c=1500;
N=10^8;
alpha=0; % absorption coefficient in dB/km

nbins=200;

x_max=5000;

%*** via the direct method
x=sqrt(rand(1,N))*x_max;

I_t_mean=100; % in dB
I_t_var=5.^2;
I_t=I_t_mean+(randn(1,N))*sqrt(I_t_var);
I_t=10.^(I_t/10);

r=sqrt(x.^2+(d-h).^2);
% in dB : Ir_dB = It_dB - 20log(r) - alpha*(r/1000). => alpha is in dB/km
% linear: Ir = It * (r^-2) * 10^(- alpha*r/(1000*10))
I_r=I_t.*(r.^(-2)); % adds spreading loss
I_r = I_r.*10.^(-alpha*r/(1000*10)); % adds absorption
%I_r=I_r+randn(1,N)*sqrt(4); % adds noise to the received intensity

P_r=sqrt(I_r).*(2*randi(2,1,length(I_r))-3);

figure
h=histogram(P_r,nbins,'normalization','pdf');
xlabel('Pressure')
ylabel('pdf')
% grid on
hold on
pd=fitdist(P_r(1:10^5).','Stable');
bins=h.BinEdges(1):(h.BinEdges(end)-h.BinEdges(1))/500:h.BinEdges(end);
f=pdf(pd,bins);
plot(bins,f);
set(gca,'yscale','log')

% minimize ||f-f_hat|| w.r.t 'x'