clc; clear; close all

% corrects in the MS sense, the divergence between the SaS PDF and
% simulated PDF by adjusting the Poisson points on the sea floor and the
% parameters of the SaS distribution

d=5;
h=10;
%c=1500;
N=10^4;
alpha=10; % absorption coefficient in dB/km

nbins=500;

x_max=5000;

%*** via the direct method
x=sqrt(rand(1,N))*x_max;
phi=2*pi*rand(1,N);
[x1,x2]=pol2cart(phi.',x.');
x_cart=[x1;x2];
%figure
%plot(x1,x2,'.','markersize',2)
%grid on
%axis equal
%pause

I_t_mean=100; % in dB
I_t_var=5.^2;
I_t=I_t_mean+(randn(1,N))*sqrt(I_t_var);
I_t=10.^(I_t/10);

r=sqrt(x.^2+(d-h).^2);
% in dB : Ir_dB = It_dB - 20log(r) - alpha*(r/1000). => alpha is in dB/km
% linear: Ir = It * (r^-2) * 10^(- alpha*r/(1000*10))
I_r=I_t.*(r.^(-2)); % adds spreading loss
I_r = I_r.*10.^(-alpha*r/(1000*10)); % adds absorption
%I_r=I_r+randn(1,N)*sqrt(4); % adds noise to the received intensity

phase=2*randi(2,1,length(I_r))-3;
P_r=sqrt(I_r).*phase; % pressure samples


edge=quantile(abs(P_r),[0, 1-10^3/N]);
edges=edge(1):(edge(2)-edge(1))/nbins:edge(2);
edges=0:600/nbins:600;

figure
histogram(abs(P_r),edges,'normalization','pdf');
xlabel('absolute pressure)')
ylabel('pdf')
% grid on
hold on
%pd=fitdist(P_r(1:min(10^5,N)).','Stable');
[a,c,~] = sstabfit(P_r);
pd = makedist('Stable','alpha',a,'beta',0,'gam',c,'delta',0);
bins = (edges(1:end-1)+edges(2:end))/2;
f = pdf(pd,bins);
plot(bins,2*f,'linewidth',2);
set(gca,'yscale','log')




% fval = cost(x_cart,d,h,alpha,I_t,edges,f);
% disp(['fval = ',num2str(fval)]);
%  
% rmse=@(x_loc)cost(x_loc,d,h,alpha,I_t,edges,f);
% [x_loc,fval] = fminunc(rmse,x_cart);
% disp(['fval = ',num2str(fval)]);
%  
% x_temp=reshape(x_cart,N,2);
% ppickingcircle(x_temp(1:min(10^4,N),:),x_max);
% x_temp=reshape(x_loc,N,2);
% ppickingcircle(x_loc(1:min(10^4,N),:),x_max)


% function J = cost(X,d,h,alpha,I_t,edges,f)
%     r = sqrt(sum(X.^2,2)+(d-h)^2);
%     Ir_dB = 10*log10(I_t.') - 20*log10(r) - alpha*(r/1000);
%     %I_r = I_t.*(r.^(-2)); % adds spreading loss
%     %I_r = I_r.*10.^(-alpha*r/(1000*10)); % adds absorption
%     I_r=10.^(Ir_dB/10);
%     p_r=sqrt(I_r); % pressure samples
%     [p_pdf,~] = histcounts(p_r,edges,'normalization','pdf');
%     J=log(sum((p_pdf-2*f).^2));
% end

function J = cost(x,d,h,alpha,I_t,edges,f)
    N=length(x)/2;
    r = sqrt(sum(reshape(x,N,2).^2,2)+(d-h)^2);
    Ir_dB = 10*log10(I_t.') - 20*log10(r) - alpha*(r/1000);
    %I_r = I_t.*(r.^(-2)); % adds spreading loss
    %I_r = I_r.*10.^(-alpha*r/(1000*10)); % adds absorption
    I_r=10.^(Ir_dB/10);
    p_r=sqrt(I_r); % pressure samples
    [p_pdf,~] = histcounts(p_r,edges,'normalization','pdf');
    J=log(sum((p_pdf-2*f).^2));
end

function ppickingcircle(X,rho)
   
    figure
    theta=0:2*pi/500:2*pi;
    [x,y]=pol2cart(theta,repmat(rho,1,length(theta)));
    plot(x,y,'-k','linewidth',2)
    hold on
    
    plot(X(:,1),X(:,2),'.','markersize',2)
    grid on
    axis equal
end