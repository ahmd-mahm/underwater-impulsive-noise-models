clc; clear; close all;

nbins=100;
samples=10^7;

xdB_min = -100;
xdB_max = 80;
x_min = db2mag(xdB_min);
x_max = db2mag(xdB_max);

xdB = xdB_min:(xdB_max-xdB_min)/nbins:xdB_max;
x = x_min:(x_max-x_min)/nbins:x_max;

pd=makedist('Stable','alpha',1.5,'beta',0,'gam',1,'delta',0);

x_rnd=abs(random(pd,samples,1));
x_rnd_dB=mag2db(x_rnd);

figure
histogram(x_rnd,x,'normalization','pdf');
grid on

figure
histogram(x_rnd_dB,xdB,'normalization','pdf');
grid on

figure
histogram(x_rnd_dB,xdB,'normalization','pdf');
set(gca,'yscale','log')
grid on
hold on
temp_x=db2mag(xdB);
g_diff_x= 10*log10(exp(1))./temp_x; % g(i)=10log10(i)
pdf_x=(pdf(pd,sqrt(temp_x))./temp_x)./g_diff_x;
plot(xdB,pdf_x)