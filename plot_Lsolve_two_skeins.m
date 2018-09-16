fonttype = 'Times';
fsize = 16;
attrib = {'Interpreter','LaTeX', ...
          'FontName',fonttype,'FontSize',fsize,'FontWeight','normal'};
lw = 2;

warning('off','Lsolve_two_skeins:Rnegative')
warning('off','Lsolve_two_skeins:dLabove1')
warning('off','Lsolve_two_skeins:dLnegative')

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
lsc = R0; tsc = 1/lambda;

% Lower bound on unraveling time (dimensional).
tlower = 1/lambda*log(Lmax/L0);

figure(1)
clf

% Solve for the two skeins case, converting to dimensionless units.
[t,L] = Lsolve_two_skeins(m,smiley,L0/R0,dLmax/R0);
% Convert back to dimensional units.
t = t*tsc; L = L*lsc;

% Upper bound on length (dimensional).
tupper = linspace(t(1),tlower,30);
Lupper = L0*exp(lambda*tupper);

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
ylabel('$L_1$ \ [$\mu$meters]',attrib{:})
set(gca,attrib{3:end})
hold off

%print -dpdf Lsolve_two_skeins
