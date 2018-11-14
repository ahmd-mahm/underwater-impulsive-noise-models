clc; clear; close all

% Impact of spreading vs aborption

r_max=100;
r=0:r_max/1000:r_max; % distance in km

sprd= -10*log10(r.^-2); % spherical spreading loss

alpha= 3; % absorption coefficient in dB/km
absp= -10*(-alpha*r/10); % absorption loss

plot(r,sprd,'linewidth',2)
hold on
plot(r,absp,'linewidth',2)
plot(r,sprd+absp,'linewidth',2) % plot total transmission loss
grid on

legend('spreading','absorption','combined loss')


ylabel('loss in dB')
xlabel('km')
grid on