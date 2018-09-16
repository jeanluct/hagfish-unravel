function [t,L,R] = Lsolve_freefree(m,P,R0,r,dLmax,unfrac)

% m, P are dimensionless.
% R0, r, dLmax scaled by L0.
% unfrac is a fraction in (0,1] denoting the unraveled thread fraction at
% which to stop.
if nargin < 6, unfrac = 1; end

L0 = 1;
tmax = inf;

f = @(t,tau) rhs(t,tau,m,P,R0,r,dLmax,unfrac);
fe = @(t,tau) events(t,tau,m,P,R0,r,dLmax,unfrac);

opts = odeset('Events',fe,'NonNegative',1,'MaxStep',1, ...
              'RelTol',1e-10,'AbsTol',1e-10);

% Solve for tau = Lmax - L.
[t,tau] = ode45(f,[0 tmax],dLmax,opts);

Lmax = L0 + dLmax;
L = Lmax - tau;
R = R0*nthroot(1 - (L-L0)/dLmax,3);

% =========================================================================
function dtau = rhs(t,tau,m,P,R0,r,dLmax,unfrac)

L0 = 1;
Lmax = L0 + dLmax;
L = Lmax - tau;

if L < 0
  warning('Lsolve_freefree:Lnegative','L=%g < 0 at t=%g.',L,t)
end

% Slendermess parameter (this is the only place we use r directly).
delta = -1/log((r/L)^2*exp(1));

if delta < 0
  warning('Lsolve_freefree:deltaneg','delta=%g < 0 at t=%g.',delta,t)
end

% R/R0
RR0 = nthroot(1 - (L-L0)/dLmax,3);

if RR0 < 0
  warning('Lsolve_freefree:Rnegative','R=%g < 0 at t=%g.',R0*RR0,t)
  dL = RR0;
  dtau = -dL;
  return
end

% The thread equation to solve for x = dL/dt.
f = @(x) x.^m + P*RR0*L*(x - .5*L)/(L + 3*R0*RR0/2/delta);

% The thread equation to solve for x = dL/dt.
dL = fsolve(f,1,optimset('Display','off','TolX',1e-15));

if dL < 0
  warning('Lsolve_freefree:dLnegative','dL/dt=%g < 0 at t=%g.',dL,t)
end

dtau = -dL;

% =========================================================================
function [value,isterm,direc] = events(t,tau,m,P,R0,r,dLmax,unfrac)

L0 = 1;
Lmax = L0 + dLmax;
L = Lmax - tau;

value(1) = unfrac - (L-L0)/dLmax;
isterm(1) = 1;
direc(1) = 0;
