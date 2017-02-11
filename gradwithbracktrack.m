%%This script performs gradient descent with backtracking line search
close all 
clear all
tic
x0=[-1.2;1];
rosenbrock=@(x1,x2)(100*(x2-x1.^2).^2+(1-x1).^2);
gradient=@(x1,x2)([-400*(x2-x1.^2)*x1-2*(1-x1);200*(x2-x1.^2)]);
gk=gradient(x0(1,1),x0(2,1));
pk=-gk;
i=1
xk=x0

%for i=1:5
while norm(gk)>=10^-3
 a_star=backtracking(xk,pk);
 a_star
 xk=xk+a_star*pk;
 gk=gradient(xk(1,1),xk(2,1));
 pk=-gk;
     i=i+1
     toc
 end