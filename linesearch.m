%%Line search algorithm that satisfies the strong wolf condition
function [astar]=linesearch(amax,x_0,pk)%%where x_k is a vector containing x1 amd x2
fprintf('Initialising line search');
c1=10^-4;%%the value was provided
c2=0.9;%%provided
a0=0;
aprev=a0;
rosenbrock=@(x1,x2)((100*(x2-x1^2).^2)+(1-x1)^2);
gradient=@(x1,x2)([-400*(x2-x1.^2)*x1-2*(1-x1);200*(x2-x1.^2)]);
a_i=amax*rand(1);
a_i 
gradphi0=(gradient(x_0(1,1),x_0(2,1)))'*pk;
gradphi_prev=gradphi0;
phi0=rosenbrock(x_0(1,1),x_0(2,1));
phi_prev=phi0;
i=1;

while (1)
     
  k=x_0+a_i*pk;
phia_i=rosenbrock(k(1,1),k(2,1));
if (phia_i>(phi0+c1*a_i*gradphi0) || (phia_i>=phi_prev && i>1))
    
    astar=zoomfun(x_0,pk,aprev,a_i);
   
    return;
end
j=x_0+a_i*pk;
gradphiai=(gradient(j(1,1),j(2,1)))'*pk;
if(abs(gradphiai)<=-c2*gradphi0)
    
    astar=a_i;
    return;
end
if(gradphiai>=0)
    astar=zoomfun(x_0,pk,a_i,aprev);
   
    return;
end

a_prev=a_i;
phi_prev=phia_i;

a_i=a_i + (amax-a_i)*rand(1);
a_i;
i=i+1;
end

    




 