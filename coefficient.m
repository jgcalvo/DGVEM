function [rho,rhs,uex,gux,guy] = coefficient(rhoC,eta)
%coefficient  Coefficient, data and exact solution of the convergence tests
%   [rho,rhs,uex,gux,guy] = coefficient(rhoC,eta) returns the coefficient 
%   rho, the right hand side f = -div(rho grad u), the exact solution and 
%   its two partial derivatives, for the five cases described in 
%   ex_convergence.m.
%
%   AUTHOR: Juan G. Calvo, 2026
uex = @(x,y) sin(pi*x).*sin(pi*y);          % default exact solution
gux = @(x,y) pi*cos(pi*x).*sin(pi*y);
guy = @(x,y) pi*sin(pi*x).*cos(pi*y);
switch rhoC
    case 0                                   % rho = 1
        rho.Fn = @(x,y) 0*x+1;   
        rho.isCont = 1;
        rhs = @(x) 2*pi*pi*sin(pi*x(:,1)).*sin(pi*x(:,2));
    case 1                                   % rho = 1 + eta x y
        rho.Fn = @(x,y) eta*x.*y+1;   
        rho.isCont = 1;
        rhs = @(x) - eta*pi*x(:,2).*cos(pi*x(:,1)).*sin(pi*x(:,2)) ...
                   - eta*pi*x(:,1).*sin(pi*x(:,1)).*cos(pi*x(:,2)) ...
          + 2*pi*pi*sin(pi*x(:,1)).*sin(pi*x(:,2)).*(eta*x(:,1).*x(:,2)+1);
    case 2                                   % rho = 2 + cos cos
        rho.Fn = @(x,y) 2+cos(10*pi*x).*cos(10*pi*y);   
        rho.isCont = 1;
        rhs = @(x) 10*pi*pi*sin(10*pi*x(:,1)).*cos(10*pi*x(:,2)).*cos(pi*x(:,1)).*sin(pi*x(:,2)) ...
                 + 10*pi*pi*cos(10*pi*x(:,1)).*sin(10*pi*x(:,2)).*sin(pi*x(:,1)).*cos(pi*x(:,2)) ...
                 + 2*pi*pi*sin(pi*x(:,1)).*sin(pi*x(:,2)).*(2+cos(10*pi*x(:,1)).*cos(10*pi*x(:,2)));
    case 3                                   % rho = 1 / eta, interface x = 1/2
        rho.Fn = @(x,y) 1+(eta-1)*(x>1/2);   
        rho.isCont = 0;
        rhs = @(x) 2*pi*pi*sin(pi*x(:,1)).*sin(pi*x(:,2)).*(1+(eta-1)*(x(:,1)>1/2));
    case 4                                   % same rho, solution with kink
        % phi(x) = x on the left; on the right a quadratic with phi
        % and rho phi' continuous at x = 1/2 and phi(1) = 0
        a   = -2-2/eta;
        phi = @(x) x.*(x<=1/2) + (a*(x-1/2).^2+(x-1/2)/eta+1/2).*(x>1/2);
        dph = @(x) 1*(x<=1/2) + (2*a*(x-1/2)+1/eta).*(x>1/2);
        rho.Fn = @(x,y) 1+(eta-1)*(x>1/2);   rho.isCont = 0;
        uex = @(x,y) phi(x).*sin(pi*y);
        gux = @(x,y) dph(x).*sin(pi*y);
        guy = @(x,y) pi*phi(x).*cos(pi*y);
        rhs = @(x) (x(:,1)<=1/2).*(pi*pi*x(:,1).*sin(pi*x(:,2))) ...
                 + (x(:,1)>1/2).*((-2*a*eta+eta*pi*pi*phi(x(:,1))).*sin(pi*x(:,2)));
end
end
