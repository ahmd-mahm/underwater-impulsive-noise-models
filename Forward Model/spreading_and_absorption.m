clc; clear; close all

% Impact of spreading vs aborption

r_max=5;
r=1/1000:r_max/1000:r_max; % distance in km

sprd= 10*log10((r*1000).^-2); % spherical spreading loss

alpha= 10; % absorption coefficient in dB/km
absp= -alpha*r; % absorption loss

% in dB : Ir_dB = It_dB - 20*log10(r) - alpha*(r/1000). => alpha is in dB/km
% linear: Ir = It * (r^-2) * 10^(- alpha*r/(1000*10))

plot(r,sprd,'linewidth',2)
hold on
plot(r,absp,'linewidth',2)
plot(r,sprd+absp,'linewidth',2) % plot total transmission loss
grid on

legend('spreading','absorption','combined loss')

title(['$\alpha=$',num2str(alpha),' dB/km'],'interpreter','latex')
ylabel('loss in dB')
xlabel('km')
grid on