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

%% *** Initial Settings ***
d=5;        % sensor depth in m
h=20;       % height of water column in m
c=1500;     % speed of sund in water in m/s
samples=10^6;   % considered time samples
fs=180000;  % samping frequency in Hz

alpha=10;   % absorption coefficient in dB/km
lambda=20;  % snaps per second
T=samples/fs; % associated time window in s

%N=round(1.1*lambda*T); % number of points (snaps)

disp(['d = ',num2str(d),' m']);
disp(['h = ',num2str(h),' m']);
disp(['T = ',num2str(T),' secs']);


%% *** ADC Noise ***

% x=loadHifDAQ3('F:\HifDAQ 2014\HIDAQ_2014\Aug_15_2014_2',0,20,1);
% x=x(:,1);
% figure
% plot(x)
% [alp,del,~]=sstabfit(x);
% y=sort(x);
% z=diff(y);
% g=z(z~=0);
% ADC_noise=20*log10(min(g));

%ADC_noise_I= 84.7; % in dB % ADC Noise
%var_noise=10.^(ADC_noise_I/10);

%% *** Evaluating x_max, tau_max, tau_min, and N ***

It_mean_dB=180; % in dB
It_var_dB=(10/3)^2;

TL_dB_max=It_mean_dB+3*sqrt(It_var_dB);
if alpha~=0
    %r=(20000*wrightOmega(it_dB/20 - Ir_dB/20 - log(20000/alpha)))/alpha;
    r_max=(20000*wrightOmega(TL_dB_max/20 - log(20000/alpha)))/alpha;
else
    r_max=10.^(TL_dB_max)/20;
end
x_max=sqrt(r_max^2-(h-d)^2);
tau_max=r_max/c;
tau_min=(h-d)/c;

Tx_interval=T-tau_min+tau_max; % transmission time window
N=round(lambda*Tx_interval); 

disp(['N = ',num2str(N),' snaps']);
disp(['x_max = ',num2str(x_max),' meters']);

%% *** PDF/CDF Transmission Loss (dB) ***
% TL_dB=it_dB-Ir_dB;
% 
% L=[10^-5,1-10^-5];
% nbins=180;
% figure
% pdfquant(TL_dB,nbins,L);
% xlabel('TL (dB) -- direct arrivals')
% 
% figure
% [f,bins]=cdfquant(TL_dB,nbins,L);
% xlabel('TL (dB) -- direct arrivals')
% 
% figure
% plot(f.Values,bins,'.','MarkerSize',16)
% hold on
% pp=pchip(f.Values,bins);
% ind_interp=f.Values(1):0.001:f.Values(end);
% pp_interp=ppval(pp,ind_interp);
% plot(ind_interp,pp_interp,'LineWidth',2)
% grid on
% ylabel('TL (dB) -- direct arrivals')
% xlabel('CDF')
% hold off
% 
% 
% TL=10.^(TL_dB/10);
% figure
% logypdfquant(TL,nbins,[10^-4,1-10^-4]);
% xlabel('TL -- direct arrivals')

%% *** Point Picking ***

x=sqrt(rand(1,N))*x_max;
phi=2*pi*rand(1,N);
x_cmp=x.*exp(1i*phi);
r=sqrt(x.^2+(h-d)^2);

t_ind_Tx=rand(1,N)*(Tx_interval)-tau_max;
t_ind_Rx= r/c+t_ind_Tx;
t_ind_Rx= t_ind_Rx(and(t_ind_Rx>=0,t_ind_Rx<=T));

ind_Rx=round(t_ind_Rx*fs);
N=length(ind_Rx);
r=r(1:N);

%% *** Generating Transmit, ADC Noise and Receive Intensities ***
It_dB=It_mean_dB+(randn(1,N))*sqrt(It_var_dB); % log-normal distribution of intensity
%it_dB=190;
Ir_dB_max=It_mean_dB+3*sqrt(It_var_dB) - 20*log10(h-d) - alpha*((h-d)/1000);

ADC_noise_lvl= (10^(Ir_dB_max/20))/(2^16); % ADC Noise with 16 bits

% in dB : Ir_dB = It_dB - 20*log10(r) - alpha*(r/1000). => alpha is in dB/km
% linear: Ir = It * (r^-2) * 10^(- alpha*r/(1000*10))

Ir_dB = It_dB - 20*log10(r) - alpha*(r/1000);
Pr = 10.^(Ir_dB/20);
Pr_ts=zeros(1,samples);
Pr_ts(ind_Rx)=Pr;

figure
stem((0:1/fs:T-1/fs),Pr_ts)
grid on
xlabel('time')
ylabel('P_r')

load('silhouette.mat');
Pr_ts=conv(Pr_ts,silh/max(silh));
figure
plot((0:length(Pr_ts)-1)/fs,Pr_ts)
grid on
xlabel('time')
ylabel('P_r (waveform)')

Ir_ts_dB=20*log10(abs(Pr_ts));
figure
plot((0:length(Pr_ts)-1)/fs,Ir_ts_dB)
grid on
xlabel('time')
ylabel('I_r (dB), waveform')

Pr_ts_ADC=Pr_ts+(rand(1,length(Pr_ts))*2-1)*ADC_noise_lvl/2;
figure
plot((0:length(Pr_ts)-1)/fs,Pr_ts_ADC)
grid on
xlabel('time')
ylabel('P_r + ADC noise (waveform)')

%% *** Histograms and PDFs ***

L=[10^-4,1-10^-4];
nbins=100;
figure
loglogpdfquant(abs(Pr_ts(Pr_ts~=0)),nbins,L);
xlabel('Pr_ts -- direct arrivals')