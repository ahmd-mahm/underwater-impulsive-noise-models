clc; close all; clear;

% Produces an orthonormal basis for the snaps recorded in 'silhouette.mat'

load('silhouette.mat','silh_mtx');
[r,c1]=size(silh_mtx);
disp(['silh_mtx : ', num2str(r),' x ',num2str(c1)]);

L=orth(silh_mtx);
[~,c2]=size(L); % reatins the same number of rows in silh_mtx, while c2=rank(silh_mtx)
disp(['orth(silh_mtx) : ', num2str(r),' x ',num2str(c2)]);

% Let 's' be a column vector (r x 1) of samples of a snap, then u=L'*s
% gives us the orthonormal weights w.r.t to the basis 'L' of 's'. The
% dimension of u is (c2 x 1). Likewise, the inequality s=L*u also holds as
% L'*s=L'*L*u=u;
%
% Vectorizing the above inequalities, we have U=L'*S and S=L*U.

S=silh_mtx;
U=L'*S;
figure(1)
p=plot(U,'.','MarkerSize',6);
grid on
xlabel('j (basis index)');
ylabel('snap coefficients')
ylim(2*[-1, 1])



for i=1:c2
    sfigure(2)
    plot(U(:,i));
    xlabel('j (basis index)');
    ylabel(['u for snap ',num2str(i)])
    ylim(2*[-1, 1])
    grid on
    
    sfigure(3)
    plot(L(:,i),'-r')
    xlabel('i (element index) ')
    ylabel(['L(:,',num2str(i),')'])
    grid on
    ylim([-1, 1])
    pause(0.1);
end