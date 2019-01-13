clc; clear; close all

% Considers a Poisson-point process (uniform distribution of events) on the
% seabed and time. The script outputs the PDF of inter-event intervals received snaps.

d=5;        % sensor depth
h=20;       % height of water column
N=10^6;     % number of points
alpha=10;   % absorption coefficient in dB/km
x_max=500;
c=1500;

SR_flag=1;

nbins =100;

rate= 15;   % observed snaps per second
T=N/rate;   % observational window of snaps

disp(['T = ',num2str(T/60),' mins']);
disp(['rate = ',num2str(rate), ' snaps per second']);

% point-picking via the direct method
x=sqrt(rand(1,N))*x_max;
phi=2*pi*rand(1,N);
x_cmp=x.*exp(1i*phi);

r_d=sqrt(x.^2+(h-d)^2);
tau_d=r_d/c; % direct arrival delay
Tx_time= rand(1,N)*T;

if SR_flag
    r_s=sqrt(x.^2+(h+d)^2);
    tau_s=r_s/c; % surface reflection delay
    Rx_time= sort([Tx_time+tau_d,Tx_time+tau_s]);
else
    Rx_time= sort(Tx_time+tau_d);
end
Rx_II= diff(Rx_time)*1000;

ppickingcircle(x_cmp(1:min(10^4,N)),x_max)

L=[10^-2,1-10^-4];

figure
pdfquant(Rx_II,nbins,L);
xlabel('II (ms)')
ylabel('pdf')

figure
logypdfquant(Rx_II,nbins,L);
xlabel('II (ms)')
ylabel('pdf')

figure
loglogpdfquant(Rx_II,nbins,L);
xlabel('II (ms)')
ylabel('pdf')


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