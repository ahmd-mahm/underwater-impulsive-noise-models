clc; clear; close all;

% Considers a Poisson-point process (uniform distribution of events) on the
% seabed and constrains the received pressure to be IID SaS. The script
% outputs the distribution of transmitted intensity.

d=5;        % sensor depth
h=20;       % height of water column
N=10^6;     % number of points
alpha=10;   % absorption coefficient in dB/km
x_max=5000;

% point-picking via the direct method
x=sqrt(rand(1,N))*x_max;
phi=2*pi*rand(1,N);
x_cmp=x.*exp(1i*phi);


%alp=1.5;    % alpha of the SaS dist.
%del=0.1;
alp=1.53;
del=2.282e+05;

pd=makedist('Stable','alpha',alp,'beta',0,'gam',del,'delta',0);
disp(['alpha = ',num2str(alp)]);
disp(['delta = ',num2str(del)]);

Pr=random(pd,1,N);
Pr_abs=abs(Pr);
Ir_dB=20*log10(Pr_abs);
Ir=10.^(Ir_dB/10);


r=sqrt(x.^2+(h-d)^2);

% in dB : Ir_dB = It_dB - 20*log10(r) - alpha*(r/1000). => alpha is in dB/km
% linear: Ir = It * (r^-2) * 10^(- alpha*r/(1000*10))

It_dB = Ir_dB + 20*log10(r)+alpha*(r/1000);
It=10.^(It_dB/10);

ppickingcircle(x_cmp(1:min(10^4,N)),x_max)

nbins = 100;
figure
[h,bins]=loglogpdfquant(Pr_abs,nbins,[10^-6,1-10^-6]);
fig_hist=gcf;
hold on

pdf_act=pdf(pd,bins);
plot(bins,2*pdf_act,'linewidth',2)
xlabel('abs(P_r)')
ylabel('pdf')

%histogram(Pr,nbins,'normalization','pdf');
%figure
%quant_It=quantile(It,[10^-4,1-10^-5]);
%bins_It= quant_It(1):(quant_It(2)- quant_It(1))/nbins:quant_It(2);
%histogram(It,bins_It);
figure
[h_It,bins_It]=loglogpdfquant(It,nbins,[10^-6,1-10^-6]);
xlabel('I_t')
ylabel('pdf')

figure
logypdfquant(It,nbins,[10^-6,1-10^-6]);
xlabel('I_t')
ylabel('pdf')

figure
dBpdfquant(It,nbins,[10^-6,1-10^-6]);
xlabel('I_t (dB)')
ylabel('pdf')

figure
dBpdfquant(Ir,nbins,[10^-6,1-10^-6]);
xlabel('I_r (dB)')
ylabel('pdf')

figure
loglogpdfquant(Ir,nbins,[10^-6,1-10^-6]);
xlabel('I_r')
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