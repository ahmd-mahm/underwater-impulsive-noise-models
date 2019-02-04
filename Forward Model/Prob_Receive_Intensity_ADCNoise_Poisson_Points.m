clc; clear; close all
% 90 dB noise / ADC res noise
% set a lambda value
% Get an impulse train for transmission times, resolve at fs
% Convolve with snap waveform.
% Get DA and SR
% Plot received intensity

% Considers a Poisson-point process (uniform distribution of events) on the
% seabed, assumes an average snap waveform, and a log-normal distribution
% of transmitted intensity. The script outputs the distribution of
% transmitted intensity.

%% *** Initial Settings ***

d=5;            % sensor depth in m
h=20;           % height of water column in m
c=1500;         % speed of sund in water in m/s
samples=10^6;   % considered time samples
fs=180000;      % samping frequency in Hz

SR=true;       % with or without surface reflections
%SR=false;

alpha=10;       % absorption coefficient in dB/km
rho=0.05;       % between 0.01 and 0.1 snaps/sec/m

T=samples/fs;   % associated time window in seconds

%N=round(1.1*lambda*T); % number of points (snaps)

disp(['d = ',num2str(d),' m']);
disp(['h = ',num2str(h),' m']);
disp(['T = ',num2str(T),' secs']);
disp(['rho = ',num2str(rho),' snaps/sec/m']);


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

%% *** Evaluating x_max, tau_max, tau_min, lambda and N ***

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


tau_min=(h-d)/c;
if SR
    r_max_sr=sqrt(x_max^2+(h+d)^2);
    tau_max=r_max_sr/c;
else
    tau_max=r_max/c;
end

Tx_interval=(T-tau_min)+tau_max; % transmission time window: [-tau_max,T-tau_min]

%*****************************
%******** LAMBDA VALUE********
%*****************************

lambda=rho*pi*(x_max^2);
%lambda=10000;

%*****************************
%*****************************

N=round(lambda*Tx_interval); % number of points (snaps)

disp(['tau_min = ',num2str(tau_min),' secs']);
disp(['tau_max = ',num2str(tau_max),' secs']);
disp(['lambda = ',num2str(lambda),' snaps/sec']);
disp(['N = ',num2str(N),' snaps in [-tau_max,T-tau_min]']);
disp(['x_max = ',num2str(x_max),' meters']);


%% *** Point Picking ***

x=sqrt(rand(1,N))*x_max;
phi=2*pi*rand(1,N);
x_cmp=x.*exp(1i*phi);
r=sqrt(x.^2+(h-d)^2);
if SR
    r_sr=sqrt(x.^2+(h+d)^2);
end

%% *** Evaluating Received Snaps' Time Indices ***

t_ind_Tx=rand(1,N)*(Tx_interval)-tau_max;
if SR
    t_ind_Rx= r/c+t_ind_Tx; 
    t_ind_Rx_sr= r_sr/c+t_ind_Tx;
    
    t_ind_logic= and(t_ind_Rx>=0,t_ind_Rx<T-1/(2*fs));      % We want a maximum length of 'samples', hence the -1/2*fs
    t_ind_sr_logic= and(t_ind_Rx_sr>=0,t_ind_Rx_sr<T-1/(2*fs));% We want a maximum length of 'samples', hence the -1/2*fs
    
%     c1=t_ind_Rx(xor(t_ind_logic,t_ind_sr_logic)); c2=t_ind_Rx_sr(xor(t_ind_logic,t_ind_sr_logic));
%     figure
%     plot(c1); hold on; plot(c2)
    
    t_ind_Rx= t_ind_Rx(t_ind_logic);
    t_ind_Rx_sr= t_ind_Rx_sr(t_ind_sr_logic);
    
    ind_Rx=round(t_ind_Rx*fs);        % resolving time indices to nearest 1/fs;
    ind_Rx_sr=round(t_ind_Rx_sr*fs);  % resolving time indices to nearest 1/fs;
    
    N=length(ind_Rx);                   % actual number of DA snaps
    N_sr=length(ind_Rx_sr);                   % actual number of DA snaps
    r=r(t_ind_logic);                   % picking DA poisson points
    r_sr=r_sr(t_ind_sr_logic);          % picking SR poisson points
else
    t_ind_Rx= r/c+t_ind_Tx;
    
    t_ind_logic= and(t_ind_Rx>=0,t_ind_Rx<T-1/(2*fs)); % We want a maximum length of 'samples', hence the -1/2*fs
    
    t_ind_Rx= t_ind_Rx(t_ind_logic);
    
    ind_Rx=round(t_ind_Rx*fs);        % resolving time indices to nearest 1/fs;
    
    N=length(ind_Rx);                   % actual number of snaps
    r=r(t_ind_logic);                   % picking 'N' poisson points
end


%% *** Generating Transmit, ADC Noise and Receive Intensities ***

% in dB : Ir_dB = It_dB - 20*log10(r) - alpha*(r/1000). => alpha is in dB/km
% linear: Ir = It * (r^-2) * 10^(- alpha*r/(1000*10))

It_dB= It_mean_dB+(randn(1,N))*sqrt(It_var_dB); % log-normal distribution of intensity
%it_dB=190;
Ir_dB_max= It_mean_dB+3*sqrt(It_var_dB) - 20*log10(h-d) - alpha*((h-d)/1000);
ADC_noise_lvl= (10^(Ir_dB_max/20))/(2^16); % ADC Noise with 16 bits

Ir_dB= It_dB - 20*log10(r) - alpha*(r/1000);
Pr= 10.^(Ir_dB/20);

Pr_ts= zeros(1,samples);
Pr_ts(ind_Rx+1)= Pr;
load('silhouette.mat','silh','silh_mtx');              % Silhouette of an average snap (DA)

if SR
    It_dB_sr= It_mean_dB+(randn(1,N_sr))*sqrt(It_var_dB); % log-normal distribution of intensity
    
    Ir_dB_sr= It_dB_sr - 20*log10(r_sr) - alpha*(r_sr/1000);
    Pr_sr= 10.^(Ir_dB_sr/20);
    
    Pr_ts_sr= zeros(1,samples);
    Pr_ts_sr(ind_Rx_sr+1)= -Pr_sr;               % negative pressure
    
    %Pr_ts_fin=conv(Pr_ts+Pr_ts_sr,silh/max(silh)); % time-series of received pressure samples
    Pr_ts_fin=snap_ts(Pr_ts+Pr_ts_sr,silh_mtx);
else
    %Pr_ts_fin= conv(Pr_ts,silh/max(silh));          % time-series of received pressure samples
    Pr_ts_fin=snap_ts(Pr_ts,silh_mtx);
end
Pr_ts_fin=Pr_ts_fin(1:samples);
Ir_ts_dB_fin= 20*log10(abs(Pr_ts_fin));

Pr_ts_ADC=Pr_ts_fin+(rand(1,samples)*2-1)*ADC_noise_lvl/2;    % adding of ADC noise


%% *** Plots ***

figure
plot((0:samples-1)/fs,Pr_ts_fin)
grid on
xlabel('time')
if SR
    ylabel('Pr (DA+SR) (waveform)')
else
    ylabel('Pr (waveform)')
end

figure
plot((0:samples-1)/fs,Ir_ts_dB_fin)
grid on
xlabel('time')
if SR
    ylabel('Ir (dB) (DA+SR), waveform')
else
    ylabel('Ir (dB), waveform')
end


figure
plot((0:length(Pr_ts_fin)-1)/fs,Pr_ts_ADC)
grid on
xlabel('time')
if SR
    ylabel('P_r + ADC noise (DA+SR) (waveform)')
else
    ylabel('P_r + ADC noise (waveform)')
end

%% *** Delay Scatter Plots ***

figure
plot(Pr_ts_ADC(1:min(end,20000)-1),Pr_ts_ADC(2:min(end,20000)),'.','markersize',1); 
grid on; 
axis equal;
xlabel('Pr_ts_ADC(i)');
ylabel('Pr_ts_ADC(i+1)')

%% *** Histograms and PDFs ***

L=[10^-5,1-10^-6];
nbins=100;
figure
loglogpdfquant(abs(Pr_ts_ADC),nbins,L);
if SR
    xlabel('Pr + ADC noise (DA+SR)')
else
    xlabel('Pr + ADC noise (DA)')
end

%% *** Point-Picking Function ***

if SR    
    %ppickingcircle(x_cmp(and(t_ind_logic,t_ind_sr_logic)),x_max,h)
    ppickingcircle(x_cmp(xor(t_ind_logic,t_ind_sr_logic)),x_max) % Plots those points, that either have a DA or a SR but not both in the received time window [0,T)
else
    %ppickingcircle(x_cmp(t_ind_logic),x_max,h)
end



%% *** Point-Picking Function ***
function ppickingcircle(x_cmp,rho,h)
if nargin==2
    figure
else
    figure(h)
end
theta=0:2*pi/500:2*pi;
[x,y]=pol2cart(theta,repmat(rho,1,length(theta)));
plot(x,y,'-k','linewidth',2)
hold on

plot(real(x_cmp),imag(x_cmp),'.','markersize',4)
grid on
axis equal
end

%% *** Snap-Waveform-Picking Function ***
function y=snap_ts(x,silh_mtx)
samples=length(x);
N=length(x~=0);

ind=(0:samples-1);
ind=ind(x~=0);

[~,c]=size(silh_mtx);

s_ind=randi(c,1,N);

for i=1:c
    temp=zeros(1,samples);
    ind_temp=ind(s_ind==i);
    temp(ind_temp+1)=x(ind_temp+1);
    
    temp=conv(temp,silh_mtx(:,i));
    temp=temp(1:samples);
    if i==1
        y=temp;
    else
        y=y+temp;
    end
end

end