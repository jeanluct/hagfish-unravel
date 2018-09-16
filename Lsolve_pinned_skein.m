function [t,L,RR0] = Lsolve_pinned_skein(m,P,r,dLmax,unfrac)

% m, P are dimensionless.
% r, dLmax scaled by L0.
% unfrac is a fraction in (0,1] denoting the unraveled thread fraction at
% which to stop.
if nargin < 5, unfrac = 1; end

L0 = 1;
tmax = inf;

f = @(t,L) rhs(t,L,m,P,r,dLmax,unfrac);
fe = @(t,L) events(t,L,m,P,r,dLmax,unfrac);

opts = odeset('Events',fe,'NonNegative',1,'RelTol',1e-6,'AbsTol',1e-6);

[t,L] = ode45(f,[0 tmax],L0,opts);

% This is R/R0, since we do not need to pass R0 to the function.
RR0 = nthroot(1 - (L-L0)/dLmax,3);

% =========================================================================
function dL = rhs(t,L,m,P,r,dLmax,unfrac)

% Slendermess parameter (this is the only place we use r directly).
delta = -1/log((r/L)^2*exp(1));

if delta < 0
  warning('Lsolve_pinned_skein:deltaneg','delta=%g < 0 at t=%g.',delta,t)
end

% The thread equation to solve for x = dL/dt.
f = @(x) (x.^m + P*L*delta*(x - 1));

dL = fsolve(f,1,optimset('Display','off','TolX',1e-15));

if dL > 1
  warning('Lsolve_pinned_skein:dLabove1','dL/dt=%g > 1 at t=%g.',dL,t)
end

if dL < 0
  warning('Lsolve_pinned_skein:dLnegative','dL/dt=%g < 0 at t=%g.',dL,t)
end

% =========================================================================
function [value,isterm,direc] = events(t,L,m,P,r,dLmax,unfrac)

L0 = 1;

value(1) = unfrac - (L-L0)/dLmax;
isterm(1) = 1;
direc(1) = 0;
