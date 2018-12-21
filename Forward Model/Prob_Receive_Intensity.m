clc; clear; close all;

% computes the PDF of the received intensity at a hydrophone receiver,
% assuming a Poisson distribution of points in time and space (sea floor).
% Considers both spreaidng and absorption losses.

d=5;
h=20;
c=1500;
N=10^7;
alpha=10; % absorption coefficient in dB/km

nbins=100;

x_max=500;

%*** via rejection sampling ***
%x=sqrt(sum(rand(2,ceil(N*4/pi)).^2))*x_max;
%x=x(x<=x_max);
%N=length(x);

%*** via the direct method
x=sqrt(rand(1,N))*x_max;



I_t_mean=100; % in dB
I_t_var=(I_t_mean*0.1/3)^2;
I_t=I_t_mean+(randn(1,N))*sqrt(I_t_var); % adds randomness to transmit intensity
%I_t=I_t_mean; % deterministic I_t

I_t=10.^(I_t/10);

%del_max=sqrt(x_max.^2+(d-h).^2)/c;
%t=(rand(1,N))*del_max;
%tim=t+r/c;

tim=rand(1,length(x));

r=sqrt(x.^2+(d-h).^2);
I_r=I_t.*(r.^(-2)); % adds spreading loss
% in dB : Ir_dB = It_dB - 20log(r) - alpha*(r/1000). => alpha is in dB/km
% linear: Ir = It * (r^-2) * 10^(- alpha*r/(1000*10))

I_r = I_r.*10.^(-alpha*r/(1000*10)); % adds absorption
%I_r=I_r+randn(1,N)*sqrt(4); % adds noise to the received intensity

p=sqrt(I_r).*(2*randi(2,1,length(I_r))-3);
pd = fitdist(p(1:min(length(p),10^5)).','Stable');

figure
histogram(x,nbins,'normalization','pdf')
xlabel('x')
ylabel('pdf')
grid on

figure
%[tim_sort,ind]=sort(tim);
%I_r=I_r(ind);
plot(tim(1:10^5),I_r(1:10^5),'.','markersize',2)
xlabel('time')
ylabel('observed intensity')
grid on


%b_max=max(I_t).*((h-d).^(-2));
if isscalar(I_t)
    b_max=I_t*((h-d).^(-2));
else
    b_max=(10^(I_t_mean/10))*((h-d).^(-2))*3;
end
%b_max=0.2;
b_min=0;
figure
histogram(I_r,(b_min:(b_max-b_min)/nbins:b_max),'normalization','pdf')
grid on
xlabel('observed intensity')
ylabel('pdf')
set(gca,'yscale','log')
grid on

hold on
%i=I_t/(x_max^2):(I_t/((h-d)^2)-I_t/(x_max^2))/1000:I_t/((h-d)^2); % theoretical PDF
%i=I_t/(x_max^2):(b_max-I_t/(x_max^2))/1000:b_max;
%f=(10^(I_t_mean/10))./((i.^2)*x_max^2);
%plot(i,f,'-k','linewidth',2)


% Analytic: Received Intensity due to deterministic I_t
i=I_t/(x_max^2):(b_max-I_t/(x_max^2))/1000:b_max;
f=(10^(I_t_mean/10))./((i.^2)*x_max^2);
plot(i,f,'-k','linewidth',2)

% Analytic: Received Intensity due to Log-Normal I_t
c1=x_max^2+(h-d)^2;
c2=(h-d)^2;
mu=I_t_mean/(10*log10(exp(1)));
sigma=sqrt(I_t_var)/(10*log10(exp(1)));

f_I=pdf('lognormal',i*c1,mu,sigma)*(c1^2/(x_max^2))+pdf('lognormal',i*c2,mu,sigma)*c2*(1-(c1/x_max^2))...
    +exp(mu+(sigma^2)/2)*sqrt(2/pi)./(2*(i.^2)*sigma*(x_max^2)) .* ( exp(-((mu+sigma^2-log(c2*i)).^2)./(2*sigma^2))...
    -exp(-((mu+sigma^2-log(c1*i)).^2)./(2*sigma^2)))...
    +exp(mu+(sigma^2)/2)./(2*(i.^2)*(x_max^2)).*(erf((mu+sigma^2-log(i*c2))/(sqrt(2)*sigma))...
    -erf((mu+sigma^2-log(i*c1))/(sqrt(2)*sigma)));
%plot(i,f_I,'-r','linewidth',2)


%F_I=cdf('lognormal',i*c1,mu,sigma)*(c1/(x_max^2))+cdf('lognormal',i*c2,mu,sigma)*(1-(c1/x_max^2))...
%    -exp(mu+(sigma^2)/2)./(2*i*(x_max^2)).*(erf((mu+sigma^2-log(i*c2))/(sqrt(2)*sigma)) - erf((mu+sigma^2-log(i*c1))/(sqrt(2)*sigma)));
