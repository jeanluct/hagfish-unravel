fonttype = 'Times';
fsize = 16;
attrib = {'Interpreter','LaTeX', ...
          'FontName',fonttype,'FontSize',fsize,'FontWeight','normal'};
lw = 2;

warning('off','Lsolve_pinned_skein:dLabove1')
warning('off','Lsolve_pinned_skein:dLnegative')

% Parameter values (dimensional)
m = 1/2;                   % force-velocity exponent
R0 = 50;                   % initial skein radius (microm)
L0 = 2*R0;                 % initial unraveled length (microm)
r = 1;                     % thread radius (microm)
eta = 1;                   % thread packing fraction in skein
dLmax = 4/3*R0^3*eta/r^2;  % maximum length added to L0 (Lmax = L0+dLmax)
Lmax = L0 + dLmax;
P = .5;

% Length and time scale.
U = 1e6;  % 1 m/s
lsc = L0; tsc = L0/U;

% Lower bound on unraveling time (dimensional).
tlower = tsc*(dLmax/lsc);

figure(1)
clf

% Solve for the pinned skein case, converting to dimensionless units.
[t,L] = Lsolve_pinned_skein(m,P,r/L0,dLmax/L0);
% Convert back to dimensional units.
t = t*tsc; L = L*lsc;

fprintf('t_unravel = %f\n',t(end));
plot(t,L,'k-','LineWidth',lw), hold on
plot([t(1) max(t(end),t(1)+tlower)],Lmax*[1 1],'r--','LineWidth',lw)
plot([t(1) t(1)+tlower],[L0 Lmax],'m:','LineWidth',lw)
axis tight
ylim([0 1.05*Lmax])
pbaspect([1 .7 1])
xlabel('$t$ \ [seconds]',attrib{:})
ylabel('$L$ \ [$\mu$meters]',attrib{:})
set(gca,attrib{3:end})
hold off

%print -dpdf Lsolve_pinned_skein
