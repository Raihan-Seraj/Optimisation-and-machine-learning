%Function for performing backtracking line search
function [a_star]=backtracking(x0,pk)
fprintf('Initialising linesearch');

abar=1;
alpha=abar;
rho=0.2;
c=10^(-4);
rosenbrock=@(x1,x2)(100*(x2-x1.^2).^2+(1-x1).^2);
gradient=@(x1,x2)([-400*(x2-x1.^2)*x1-2*(1-x1);200*(x2-x1.^2)]);
l=x0+alpha*pk
fkalphapk=rosenbrock(l(1,1),l(2,1));
fk=rosenbrock(x0(1,1),x0(2,1));
while(fkalphapk>fk+c*alpha*(gradient(x0(1,1),x0(2,1)))'*pk)

    alpha=rho*alpha;
    l=x0+alpha*pk;
    fkalphapk=rosenbrock(l(1,1),l(2,1));
    
end
a_star=alpha;
    
