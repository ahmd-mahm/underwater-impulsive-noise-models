clc; clear; close all;

% computes the PDF of the received intensity at a hydrophone receiver,
% assuming a Poisson distribution of points in time and space (sea floor).
% Considers both spreaidng and absorption losses.

d=10;
h=20;
c=1500;
N=10^8;
alpha=5; % absorption coefficient in dB/km

nbins=100;

x_max=500;

%*** via rejection sampling ***
%x=sqrt(sum(rand(2,ceil(N*4/pi)).^2))*x_max;
%x=x(x<=x_max);
%N=length(x);

%*** via the direct method
x=sqrt(rand(1,N))*x_max;



%I_t=10000*(rand(1,N).^2);
I_t=10000+(2*rand(1,N)-1)*100000*0.1;

%del_max=sqrt(x_max.^2+(d-h).^2)/c;
%t=(rand(1,N))*del_max;
%tim=t+r/c;



tim=rand(1,length(x));

r=sqrt(x.^2+(d-h).^2);
I_r=I_t.*(r.^(-2)); % in dB : Ir_dB = It_dB - 20log(r) - alpha*(r/1000). => alpha is in dB/km 
                    % linear: Ir = It * (r^-2) * 10^(- alpha*r/(1000*10))

I_r = I_r.*10.^(-alpha*r/(1000*10));
I_r=I_r+randn(1,N)*sqrt(4);

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


b_max=max(I_t).*((h-d).^(-2));
%b_max=0.2;
b_min=0;
figure
histogram(I_r,(b_min:(b_max-b_min)/nbins:b_max),'normalization','pdf')
grid on
xlabel('observed intensity')
ylabel('pdf')
set(gca,'yscale','log')
grid on