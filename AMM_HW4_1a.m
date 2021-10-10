%given parameters
c0 = 343;
rho0 = 1.21;
ln = .011;
lc = .0535;
Dn = .012;
Dc = .05;
Swg = .01;
Sn = pi*(Dn/2)^2;
Vc = pi*(Dc/2)^2*lc;
lp = 4*Dn/(3*pi);
Z0 = c0*rho0/Swg;    

%variable frequency
f = 0:1:100;
omega = (2*pi).*f;
k = omega./c0;

%real and imaginary parts of Z_HR
% %...Blackstock's Equation for the HR
% ZHR = rho0*c0.*k.^2./(2*pi) + 1i*(omega.*rho0).*(lp/Sn-(c0^2./(Vc.*omega.^2)));

%...Homeowork Equation for the Hr
ZHR = rho0*c0.*k.^2./(2*pi) + 1i*(omega.*rho0./Sn).*(ln+2*lp-(c0*Sn./(Vc.*omega.^2)));
f0 = Sn/(2*pi)*sqrt(c0/(Vc*(ln+2*lp)));
R = -Z0./(Z0+2.*ZHR);
T = 2.*ZHR./(Z0+2.*ZHR);

reR = real(R);
imR = imag(R);
reT = real(T);
imT = imag(T);
magR = abs(R);
magT = abs(T);

%plotting
tiledlayout(1,3)

nexttile
plot(f,reR,f,imR)
legend('Re $R$','Im $R$','Interpreter','latex')
hold on
grid

nexttile
plot(f,reT,f,imT)
legend('Re $T$','Im $T$','Interpreter','latex')
grid

nexttile
plot(f,magR,f,magT)
legend('$|R|$','$|T|$','Interpreter','latex')
grid
hold off

clear
