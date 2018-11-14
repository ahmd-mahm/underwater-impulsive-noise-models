clc; clear; close all;

% computes the PDF of the received intensity at a hydrophone receiver,
% assuming a Poisson distribution of points in time and space (sea floor).
% Considers both spreaidng and absorption losses.

d=10;
h=20;
c=1500;
N=10^8;
alpha=10;

x_max=5000;

x=sqrt(sum(rand(2,ceil(N*4/pi)).^2))*x_max;
x=x(x<=x_max);
N=length(x);

%I_t=10000*(rand(1,N).^2);
I_t=10000+(2*rand(1,N)-1)*100000*0.05*0;

%del_max=sqrt(x_max.^2+(d-h).^2)/c;
%t=(rand(1,N))*del_max;
%tim=t+r/c;



tim=rand(1,length(x));

r=sqrt(x.^2+(d-h).^2);
I_r=I_t.*(r.^(-2)); % in dB : Ir_dB = It_dB - 20log(r) - alpha*(r/1000). => alpha is in dB/km 
                    % linear: Ir = It * (r^-2) * 10^(- alpha*r/(1000*10))

I_r = I_r.*10.^(-alpha*r/(1000*10));


figure
%[tim_sort,ind]=sort(tim);
%I_r=I_r(ind);
plot(tim(1:10^5),I_r(1:10^5),'.','markersize',2)
xlabel('time')
ylabel('observed intensity')
grid on


b_max=max(I_t).*((h-d).^(-2));
%b_max=0.2;
b_min=0;
nbins=100;
figure
histogram(I_r,(b_min:(b_max-b_min)/nbins:b_max),'normalization','pdf')
hold on
nbins=200;
histogram(I_r,(b_min:(b_max-b_min)/nbins:b_max),'normalization','pdf')
xlabel('observed intensity')
ylabel('pdf')
set(gca,'yscale','log')
grid on