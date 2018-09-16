fonttype = 'Times';
fsize = 16;
attrib = {'Interpreter','LaTeX', ...
          'FontName',fonttype,'FontSize',fsize,'FontWeight','normal'};
lw = 2;

warning('off','Lsolve_freefree:Rnegative')
warning('off','Lsolve_freefree:dLnegative')

% Parameter values (dimensional)
m = 1/2;                   % force-velocity exponent
R0 = 50;                   % initial skein radius (microm)
L0 = 2*R0;                 % initial unraveled length (microm)
r = 1;                     % thread radius (microm)
eta = 1;                   % thread packing fraction in skein
dLmax = 4/3*R0^3*eta/r^2;  % maximum length added to L0 (Lmax = L0+dLmax)
Lmax = L0 + dLmax;
smiley = 10;

% Length and time scale.
lambda = 10;  % 10 s
lsc = L0; tsc = 1/lambda;

% Lower bound on unraveling time (dimensional).
tlower = 2/lambda*log(Lmax/L0);

figure(1)
clf

% Solve for the pinned skein case, converting to dimensionless units.
[t,L] = Lsolve_freefree(m,smiley,R0/L0,r/L0,dLmax/L0);
% Convert back to dimensional units.
t = t*tsc; L = L*lsc;

% Upper bound on length (dimensional).
tupper = linspace(t(1),tlower,30);
Lupper = L0*exp(.5*lambda*tupper);

pfun = @plot;  % set pfun = @plot or pfun = @loglog to change plot type.
if isequal(pfun,@loglog), t = t(2:end); L = L(2:end); end

fprintf('t_unravel = %f\n',t(end));
pfun(t,L,'k-','LineWidth',lw), hold on
pfun([t(1) max(t(end),t(1)+tlower)],Lmax*[1 1],'r--','LineWidth',lw)
pfun(tupper,Lupper,'m:','LineWidth',lw)
axis tight
ylim([0 1.05*Lmax])
pbaspect([1 .7 1])
xlabel('$t$ \ [seconds]',attrib{:})
ylabel('$L$ \ [$\mu$meters]',attrib{:})
set(gca,attrib{3:end})
hold off

%print -dpdf Lsolve_freefree

figure(2)
ell = (Lmax - L)/lsc; ell = ell(1:end-1);
T = (t(end) - t)/tsc; T = T(1:end-1);
mm = (3*m-1)/(3*m);
C = (.5*smiley*(Lmax/L0)/(Lmax/L0 - 1)^(1/3))^(1/m);
ellT = (mm*C*T).^(1/mm);
loglog(T,ell,'k.-')
hold on
loglog(T,ellT,'b--')
hold off
pbaspect([1 .7 1])
xlabel('$T$',attrib{:})
ylabel('$\ell$',attrib{:})
set(gca,attrib{3:end})
axis tight
ylim([min(ell) max(ell)])
