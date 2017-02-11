close all
clear all
clc
%%In this script we implement the BFGS algortihm on rosenbrock functions


rosenbrock=@(x1,x2)((100*(x2-x1^2).^2)+(1-x1)^2);
gradient=@(x1,x2)([-400*(x2-x1.^2)*x1-2*(1-x1);200*(x2-x1.^2)]);
for u=1:2
 if(u==1)
 x_0=[1.2;1.2];
 else
 x_0=[-1.2;1];
 end
x_prev=x_0;
x_k=x_0;
H0=eye(2)%% initialising H_0

g0=gradient(x_0(1,1),x_0(2,1));
gk=g0
tolerance=10^-3;
Hk=H0;
pk=-H0*g0
g_prev=g0;

i=1;

fprintf('Initialising BGGS algorithm');
while(norm(gk)>=tolerance)
    
   alpha=linesearch(10,x_k,pk);
   
   x_k=x_prev+alpha*pk;
   gk=gradient(x_k(1,1),x_k(2,1));
   s_k=x_k-x_prev;
   y_k=gk-g_prev;
   rho_k= 1/((y_k)'*s_k);
   Hk=(eye(2)-rho_k*s_k*y_k')*Hk*(eye(2)-rho_k*y_k*s_k')+rho_k*s_k*s_k';
    pk=-Hk*gk;
   g_prev=gk;
   
   x_prev=x_k;
   
   cost(i,u)=rosenbrock(x_k(1,1),x_k(2,1));
   x1vec(i,u)=x_k(1,1);
   x2vec(i,u)=x_k(2,1);
   i=i+1;
   
end
end
 
% fprintf('\nvalue of the function of x0=[1.2,1.2]\n %d',cost(end,1));
% fprintf('\nvalue of the function of x0=[-1.2,1]\n %d',cost(end,1));
%=======================================================================================================
%plotting the contours
%%==============================================================================================================
%Contour Plots

       
 f = @(x) 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
xrange = linspace( -1.5,1.5,200);
yrange = linspace( -0.5,1.5,200);
[A,B] = meshgrid( xrange, yrange );
vecA = A(:); vecB = B(:);
 for i = 1:length(vecA)
z(i,1) = f( [vecA(i); vecB(i)] );
 end
Z = reshape(z, 200, 200);

%Plot steps
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2-200 scrsz(3)/1.2 scrsz(4)/2+100])

%set linewidth and fontsize
 set(0,'defaultLineLineWidth',1.5);
set(0,'defaultAxesFontSize',16)

subplot(1,2,1)
 contourf(A,B,Z,20,'linestyle','None');
axis square
hold on
plot( x1vec(:,1)', x2vec(:,1), '-wo', 'linewidth',2)
 xlabel('x_1')
ylabel('x_2')
title(sprintf('BFGS, x_0 =[1.2,1.2]'));

 subplot(1,2,2)
 contourf(A,B,Z,20,'linestyle','None');
axis square
hold on
plot( x1vec(:,2)', x2vec(:,2), '-wo', 'linewidth',2)
 xlabel('x_1')
ylabel('x_2')
title(sprintf('BFGS, x_0 =[-1.2,1]'));

 
    
    
  