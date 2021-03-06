function [t,L,R] = Lsolve_two_skeins(m,P,L0,dLmax,unfrac)

% m, P are dimensionless.
% L0, dLmax scaled by R0.
% unfrac is a fraction in (0,1] denoting the unraveled thread fraction at
% which to stop.
if nargin < 5, unfrac = 1; end

R0 = 1;
tmax = inf;

f = @(t,tau) rhs(t,tau,m,P,L0,dLmax,unfrac);
fe = @(t,tau) events(t,tau,m,P,L0,dLmax,unfrac);

opts = odeset('Events',fe,'RelTol',1e-8,'AbsTol',1e-8,'NonNegative',1);

% Solve for tau = Lmax - L.
[t,tau] = ode45(f,[0 tmax],dLmax,opts);

Lmax = L0 + dLmax;
L = Lmax - tau;
R = R0*nthroot(1 - (L-L0)/dLmax,3);

% =========================================================================
function dtau = rhs(t,tau,m,P,L0,dLmax,unfrac)

Lmax = L0 + dLmax;
L = Lmax - tau;
R = nthroot(1 - (L-L0)/dLmax,3);

if R < 0
  warning('Lsolve_two_skeins:Rnegative','R=%g < 0 at t=%g.',R,t)
end

f = @(x) (x.^m + P*R*(x - L));

% The thread equation to solve for x = dL/dt.
dL = fsolve(f,1,optimset('Display','off','TolX',1e-15));

if dL > 1
  warning('Lsolve_two_skeins:dLabove1','dL/dt=%g > 1 at t=%g.',dL,t)
end

if dL < 0
  warning('Lsolve_two_skeins:dLnegative','dL/dt=%g < 0 at t=%g.',dL,t)
end

dtau = -dL;

% =========================================================================
function [value,isterm,direc] = events(t,tau,m,P,L0,dLmax,unfrac)

Lmax = L0 + dLmax;
L = Lmax - tau;

value(1) = unfrac - (L-L0)/dLmax;
isterm(1) = 1;
direc(1) = 0;
