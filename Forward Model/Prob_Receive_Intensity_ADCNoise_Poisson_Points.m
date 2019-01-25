clc; clear; close all

% 90 dB noise / ADC res noise
% set a lambda value
% Get an impulse train for transmission times, resolve at 180 kHz 
% Convolve with snap waveform.
% Get DA and SR
% Plot received intensity



% Considers a Poisson-point process (uniform distribution of events) on the
% seabed, assumes an average snap waveform, and a log-normal distribution
% of transmitted intensity. The script outputs the distribution of
% transmitted intensity.

d=5;        % sensor depth
h=20;       % height of water column
N=10^7;     % number of points
alpha=10;   % absorption coefficient in dB/km
x_max=200;

lambda=20; % snaps per second

% x=loadHifDAQ3('F:\HifDAQ 2014\HIDAQ_2014\Aug_15_2014_2',0,20,1);
% x=x(:,1);
% figure
% plot(x)
% [alp,del,~]=sstabfit(x);
% y=sort(x);
% z=diff(y);
% g=z(z~=0);
% ADC_noise=20*log10(min(g));

ADC_noise_I= 84.7; % in dB

alp=1.53;
del=2.282e+05;
pd=makedist('Stable','alpha',alp,'beta',0,'gam',del,'delta',0);
disp(['alpha = ',num2str(alp)]);
disp(['delta = ',num2str(del)]);

Pr=random(pd,1,N);
Pr_abs=abs(Pr);
Ir_dB=20*log10(Pr_abs);

x_max=


% in dB : Ir_dB = It_dB - 20*log10(r) - alpha*(r/1000). => alpha is in dB/km
% linear: Ir = It * (r^-2) * 10^(- alpha*r/(1000*10))

It_mean_dB=250; % in dB
It_var_dB=(10/3)^2;
it_dB=It_mean_dB+(randn(1,N))*sqrt(It_var_dB); % log-normal distribution of intensity
%it_dB=190;

TL_dB=it_dB-Ir_dB;

L=[10^-5,1-10^-5];
nbins=180;
figure
pdfquant(TL_dB,nbins,L);
xlabel('TL (dB) -- direct arrivals')

figure
[f,bins]=cdfquant(TL_dB,nbins,L);
xlabel('TL (dB) -- direct arrivals')

figure
plot(f.Values,bins,'.','MarkerSize',16)
hold on
pp=pchip(f.Values,bins);
ind_interp=f.Values(1):0.001:f.Values(end);
pp_interp=ppval(pp,ind_interp);
plot(ind_interp,pp_interp,'LineWidth',2)
grid on
ylabel('TL (dB) -- direct arrivals')
xlabel('CDF')
hold off


TL=10.^(TL_dB/10);
figure
logypdfquant(TL,nbins,[10^-4,1-10^-4]);
xlabel('TL -- direct arrivals')

points=10^4;
U=rand(1,2*points);
U=U(and(U<=f.Values(end),U>=f.Values(1)));
U=U(1:points);
TL_dB_rand=ppval(pp,U);

%syms Ir_dB It_dB r alpha
%eqn = Ir_dB == It_dB - 20*log(r) - alpha*(r/1000);
%r= solve(eqn,r)

if alpha~=0
    %r=(20000*wrightOmega(it_dB/20 - Ir_dB/20 - log(20000/alpha)))/alpha;
    r=(20000*wrightOmega(TL_dB_rand/20 - log(20000/alpha)))/alpha;
else
    r=10.^(TL_dB_rand)/20;
end
r=r(r>h-d);
points=length(r);

x=sqrt(r.^2-(h-d)^2);
phi=2*pi*rand(1,points);
x_cmp=x.*exp(1i*phi);
% 
ppickingcircle(x_cmp(1:min(10^4,points)),x_max)
title(['Total Points = ', num2str(points)])

figure
[~,bins]=loglogpdfquant(Pr_abs,nbins,L);
fig_hist=gcf;
hold on
pdf_act=pdf(pd,bins);
plot(bins,2*pdf_act,'linewidth',2)
xlabel('abs(P_r) -- direct arrivals')


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