%given parameters
c0 = 343;
rho0 = 1.21;
ln = .011;
lc1 = .0535;
lc2 = .06;
Dn = .012;
Dc = .05;
Swg = .01;
Sn = pi*(Dn/2)^2;
Vc1 = pi*(Dc/2)^2*lc1;
Vc2 = pi*(Dc/2)^2*lc2;
lp = 8*(Dn/2)/(3*pi);
L = .056;
Z0 = c0*rho0/Swg;

%variable frequency
f = 0:.01:600;
omega = (2*pi).*f;
k = omega./c0;

% Values for the first Helmholtz Resonator

%...Blackstocks equations
% ZHR1 = rho0*c0.*k.^2./(2*pi) + 1i*(omega.*rho0).*(lp/Sn-(c0^2./(Vc1.*omega.^2)));
% ZHR2 = rho0*c0.*k.^2./(2*pi) + 1i*(omega.*rho0).*(lp/Sn-(c0^2./(Vc2.*omega.^2)));

% %...Homework equations
ZHR1 = rho0*c0.*k.^2./(2*pi) + 1i*(omega.*rho0./Sn).*(ln+2*lp-(c0*Sn./(Vc1.*omega.^2)));
ZHR2 = rho0*c0.*k.^2./(2*pi) + 1i*(omega.*rho0./Sn).*(ln+2*lp-(c0*Sn./(Vc2.*omega.^2)));

%values of R and T

R1 = -Z0./(Z0+2.*ZHR1);
T1 = 2.*ZHR1./(Z0+2.*ZHR1);
R2 = -Z0./(Z0+2.*ZHR2);
T2 = 2.*ZHR2./(Z0+2.*ZHR2);

S11 = zeros(1,length(f));
S12 = zeros(1,length(f));
S21 = zeros(1,length(f));
S22 = zeros(1,length(f));

%constituent M matrices
for x = 1:length(f)
    M1 = [T1(x)-R1(x)^2/T1(x) R1(x)/T1(x);-R1(x)/T1(x) 1/T1(x)];
    M2 = [T2(x)-R2(x)^2/T2(x) R2(x)/T2(x);-R2(x)/T2(x) 1/T2(x)];
    M3 = [exp(-1i*k(x)*L) 0;0 exp(1i*k(x)*L)];
    
    Mtot = M1*M3*M2;
    
    S11(x) = -1*Mtot(1,1)*Mtot(2,1)/(1+Mtot(1,2)*Mtot(2,1));
    S12(x) = 1./Mtot(2,2);
    S21(x) = Mtot(1,1)/(1+Mtot(1,2)*Mtot(2,1));
    S22(x) = Mtot(1,2)/Mtot(2,2);
end

%convert from M to S matrix

tiledlayout(2,3)
nexttile
plot(f,real(S11),f,real(S12),f,real(S21),f,real(S22))
legend('$Re S_{11}$','$Re S_{21}$','$Re S_{12}$','$Re S_{22}$','Interpreter','latex')
xlim([0 100])
grid

nexttile
plot(f,imag(S11),f,imag(S12),f,imag(S21),f,imag(S22))
legend('$Im S_{11}$','$Im S_{21}$','$Im S_{12}$','$Im S_{22}$','Interpreter','latex')
xlim([0 100])
grid

nexttile
plot(f,abs(S11),f,abs(S12),f,abs(S21),f,abs(S22))
legend('$|S_{11}|$','$|S_{21}|$','$|S_{12}|$','$|S_{22}|$','Interpreter','latex')
xlim([0 100])
grid

nexttile
plot(f,real(S11),f,real(S12),f,real(S21),f,real(S22))
legend('$Re S_{11}$','$Re S_{21}$','$Re S_{12}$','$Re S_{22}$','Interpreter','latex')
xlim([200 600])
grid

nexttile
plot(f,imag(S11),f,imag(S12),f,imag(S21),f,imag(S22))
legend('$Im S_{11}$','$Im S_{21}$','$Im S_{12}$','$Im S_{22}$','Interpreter','latex')
xlim([200 600])
grid

nexttile
plot(f,abs(S11),f,abs(S12),f,abs(S21),f,abs(S22))
legend('$|S_{11}|$','$|S_{21}|$','$|S_{12}|$','$|S_{22}|$','Interpreter','latex')
xlim([200 600])
grid

hold off
clear
