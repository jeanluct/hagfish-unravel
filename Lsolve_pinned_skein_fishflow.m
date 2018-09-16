function [t,L,RR0] = Lsolve_pinned_skein_fishflow(m,P,r,dLmax,unfrac)

% m, P are dimensionless.
% r, dLmax scaled by L0.
% unfrac is a fraction in (0,1] denoting the unraveled thread fraction at
% which to stop.
if nargin < 5, unfrac = 1; end

L0 = 1;
tmax = inf;

f = @(t,L) rhs(t,L,m,P,r,dLmax,unfrac);
fe = @(t,L) events(t,L,m,P,r,dLmax,unfrac);

opts = odeset('Events',fe,'NonNegative',1,'RelTol',1e-10,'AbsTol',1e-10);

[t,L] = ode45(f,[0 tmax],L0,opts);

% This is R/R0, since we do not need to pass R0 to the function.
RR0 = nthroot(1 - (L-L0)/dLmax,3);

% =========================================================================
function dL = rhs(t,L,m,P,r,dLmax,unfrac)

% Slendermess parameter (this is the only place we use r directly).
delta = -1/log((r/L)^2*exp(1));

gdia = 1e5;  % gape size of predator (in microns)
xfish = 3e4; % distance of gape from the pinned skein location (in microns)

% L is converted here in dimensional units (microns), only for the purpose
% of calculating velocity field
check = 100*L;

%
% Suction velocity (Eq. 3.1 in Supplementary Information)
%
velocity = @(Lx) ...
    0.098*((100*Lx-xfish)/(gdia)).^4 + ...
    0.700*((100*Lx-xfish)/(gdia)).^3 + ...
    1.860*((100*Lx-xfish)/(gdia)).^2 + ...
    2.190*((100*Lx-xfish)/(gdia)) + 1;
% Check if the unraveled length is less than the distance between fish and
% gape.
if check < xfish
  avgvelocity = (1/L)*integral(velocity,0,L); 
else
  avgvelocity = (1/L)*(integral(velocity,0,(xfish/100)) + 1*(L-(xfish/100)));
end

if delta < 0
  warning('Lsolve_pinned_skein:deltaneg','delta=%g < 0 at t=%g.',delta,t)
end

% The thread equation to solve for x = dL/dt.
f = @(x) (x.^m + P*L*delta*(x - avgvelocity));

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
