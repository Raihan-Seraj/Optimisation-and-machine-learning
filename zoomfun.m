function astar = zoomfun(x0,pk,alphal,alphah)

c1 = 10^-4;
c2 = 0.9;
rosenbrock=@(x1,x2)((100*(x2-x1^2).^2)+(1-x1)^2);
gradient=@(x1,x2)([-400*(x2-x1.^2)*x1-2*(1-x1);200*(x2-x1.^2)]);
fx0=rosenbrock(x0(1,1),x0(2,1));
gx0=(gradient(x0(1,1),x0(2,1)))'*pk;


while(1)
   aj = 1/2*(alphal+alphah);
   k=x0+aj*pk;
   phiaj=rosenbrock(k(1,1),k(2,1));
   phi0=rosenbrock(x0(1,1),x0(2,1));
   gradphi0=(gradient(x0(1,1),x0(2,1)))'*pk;
   
   z=x0+alphal*pk;
   phialow=rosenbrock(z(1,1),z(2,1));
   
   
   %fxl = feval(f,xl,d);
   if ((phiaj > phi0 + c1*aj*gradphi0) || (phiaj >= phialow)),
       
      
      alphah = aj;
      
      
      
   else
        
        k=x0+aj*pk;
        gradphiaj=(gradient(k(1,1),k(2,1)))'*pk;
      if abs(gradphiaj) <= -c2*gradphi0
          
        astar = aj;
        return;
      end
      if gradphiaj*(alphah-alphal) >= 0
        alphah = alphal;
      end
      alphal = aj;
   end
end 