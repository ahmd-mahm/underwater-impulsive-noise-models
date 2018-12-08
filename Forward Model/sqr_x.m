clc; clear; close all;

nbins=100;
samples=10^7;

ydB_min = -100;
ydB_max = 80;
x_min = db2mag(ydB_min);
x_max = db2mag(ydB_max);

ydB = ydB_min:(ydB_max-ydB_min)/nbins:ydB_max;
x = x_min:(x_max-x_min)/nbins:x_max;

pd=makedist('Stable','alpha',1.5,'beta',0,'gam',1,'delta',0);

x_rnd=abs(random(pd,samples,1));
y_rnd_dB=mag2db(x_rnd);

figure
y_rnd=x_rnd.^2;
%y=x_min^2:(x_max^2-x_min^2)/nbins:x_max^2;
y_int=quantile(y_rnd,[10^-4, 1-10^-4 ]);
y= y_int(1):(y_int(2)-y_int(1))/nbins:y_int(2);
%y=0:0.01/nbins:0.01;
histogram(y_rnd,y,'normalization','pdf');
grid on
hold on
y_bins=(y(1:end-1)+y(2:end))/2;
pdf_y= pdf(pd,sqrt(y_bins))./sqrt(y_bins); % y=x^2;
plot(y_bins,pdf_y,'linewidth',2)
set(gca,'yscale','log')

% figure
% histogram(y_rnd_dB,ydB,'normalization','pdf');
% grid on
% 
figure
histogram(y_rnd_dB,ydB,'normalization','pdf');
set(gca,'yscale','log')
grid on
hold on
y=db2pow(ydB);
pdf_y= pdf(pd,sqrt(y))./sqrt(y); % y=x^2
g_diff_y= 10*log10(exp(1))./y; % g(y)=10log10(y)
pdf_y_dB=pdf_y./g_diff_y;
plot(ydB,pdf_y_dB,'linewidth',2)
G=gradient(log10(pdf_y_dB)*10,ydB);
G(1:2)