%%This script performs newtons method using backtracking line search 
close all
clear all
clc
%newton's method 
tic
rosenbrock=@(x1,x2)(100*(x2-x1.^2).^2+(1-x1).^2);
gradient=@(x1,x2)([-400*(x2-x1.^2)*x1-2*(1-x1);200*(x2-x1.^2)]);
 hessian=@(x1,x2)([-400*(x2-x1.^2)+800*x1.^2+2 -400*x1;-400*x1 200]);
 x0=[-1.2;1];
 bk=hessian(x0(1),x0(2));
 gradk=gradient(x0(1),x0(2));
 pk=-1*inv(bk)*gradk;
xk=x0;
i=0
 while norm(gradk)>=10^-3
           alpha=backtracking(xk,pk)
            
            xk=xk+alpha*pk;
            bk=hessian(xk(1,1),xk(2,1));
            gradk=gradient(xk(1,1),xk(2,1));
            pk=-1*inv(bk)*gradk;
            
            
      i=i+1      
      toc
     
  end