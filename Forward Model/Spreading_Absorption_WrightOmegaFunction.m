clc; clear; close all

% Ascertains the valdity of the Wright-Omega solution for deducing the
% range in the presence of geometric spreading and absorption, for given
% transmit and receive intensities, and absorption coefficient.


%syms Ir_dB It_dB r alpha
%eqn = Ir_dB == It_dB - 20*c*log(r) - alpha*(r/1000);
%r= solve(eqn,r)
%where c= 1/log(10)
%r = (1000*c*wrightOmega(- log((1000*c)/alpha) - (1000*Ir_dB - 1000*It_dB)/(1000*c)))/alpha

% in dB : Ir_dB = It_dB - 20*log10(r) - alpha*(r/1000). => alpha is in dB/km
% linear: Ir = It * (r^-2) * 10^(- alpha*r/(1000*10))

It_dB=180;
Ir_dB=75;

r_max=1000;
amin=-20; amax=10;
a_ind=10.^(amin:(amax-amin)/100:amax);

r=zeros(1,length(a_ind));
flag=zeros(1,length(a_ind));
sprd_only_range= 10^((It_dB - Ir_dB)./20);

i=1;
for a=a_ind
 %r(i)=(20000*wrightOmega(It_dB/20 - Ir_dB/20 - log(20000/a_ind(i))))/a_ind(i);
 r(i)=(1000*20/log(10)*wrightOmega(- log((1000*20/log(10))/a_ind(i)) - (1000*Ir_dB - 1000*It_dB)/(1000*20/log(10))))/a_ind(i);
 flag(i)= (It_dB - Ir_dB)./(20*log10(r(i)) + a*(r(i)/1000));
 i=i+1;
end

figure
yyaxis left
semilogx(a_ind,r,'linewidth',2)
hold on
semilogx(a_ind,repmat(sprd_only_range,1,length(a_ind)),'--k')
ylabel('r_{max}')
xlabel('alpha (dB/km)')
grid on

yyaxis right
plot(a_ind,flag)
ylim([-1, 1]*0.5+1)
ylabel('Sanity Check')

legend('r_{max} (with spreading + absorption)','r_{max} (with spreading)','Validity of Wright-Omega Solution (==1)','Location','southwest')