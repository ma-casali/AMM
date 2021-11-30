%Format of ForwardPress and BackwardPress:
%f,rep1_s,rep2_s,rep1_b,rep2_b,imp1_s,imp2_s,imp1_b,imp2_b
%ForwardPress:
%columns 4,6,9,11 can be expelled
%BackwardPress:
%columns 3,7,9 can be expelled

c0 = 343; %m/s
rho0 = 1.15; %kg/m^3
Z0 = c0*rho0; %rayls
d12 = .00224; %m
d2s = .00183; %m
d3s = .00183; %m
d34 = .00224; %m
L = 0.006475; %m (length of inclusion)

%Curate Matrices of Data
f_press = forwardpress;
b_press = backwardpress;

f_press(:,[4 6 10]) = [];
b_press(:,[3 4 8]) = [];

freq = f_press(:,1);
omega = 2*pi.*freq;
k0 = omega./c0;

for x = 1:length(freq)
    R(x,1) = (f_press(x,2)+1i*f_press(x,6))/(f_press(x,4)+1i*f_press(x,8));
    T(x,1) = abs(f_press(x,3)+1i*f_press(x,7))*exp(-1i*(d3s+d34)*k0(x)/2);
    
    R_B(x,1) = (b_press(x,3)+1i*b_press(x,7))/(b_press(x,5)+1i*b_press(x,9));
    T_B(x,1) = abs(b_press(x,2)+1i*b_press(x,6))*exp(-1i*(d3s+d34)*k0(x)/2);
end

%Plotting Magnitude and Phase of the Reflection and Transmission Coefficients
f1 = figure();
figure(f1);
tiledlayout(1,2)
nexttile
plot(freq,abs(R),'s',freq,abs(R_B),'x')
hold on
plot(freq,abs(T),'s',freq,abs(T_B),'x')
hold off
legend("R","R_B","T","T_B")
xlabel("Frequency (Hz)")
title("Magnitude")
grid

nexttile
plot(freq,imag(R),'s',freq,imag(R_B),'x',freq,imag(T),'s',freq,imag(T_B),'x')
legend("R","R_B","T","T_B")
xlabel("Frequency (Hz)")
title("Phase")
grid

%Effective Material Properties -------------------------------------------
%Z,k and W
m = floor(L./(freq./c0));
r = sqrt((1-R.*R_B+T.^2).^2 - T.^2);
x = (1-R.*R_B+T.^2+r)./(2.*T);

Z = (Z0.*r)./((1-R).*(1-R_B)-T.^2);
k = 1/L.*(1i.*log(x)+2*pi.*m);
W = (R_B - R)./(1i.*r);

kappa = Z.*omega./k;
rho = Z.*k./omega;
psi = W.*Z./omega;

%plotting effective properties
f2 = figure();
figure(f2);
tiledlayout(2,2)

nexttile
plot(freq,real(rho),'s',freq,imag(rho),'^')
grid
legend("Real","Imag")
xlabel("Frequency (Hz)")
ylabel("$\rho$ ($kg/m^3$)","Interpreter","Latex")

nexttile
plot(freq,real(kappa)/1000,'s',freq,imag(kappa)/1000,'^')
grid
legend("Real","Imag")
xlabel("Frequency (Hz)")
ylabel("$\kappa$ (kPa)","Interpreter","Latex")

nexttile
plot(freq,real(psi),'s',freq,imag(psi),'^')
grid
legend("Real","Imag")
xlabel("Frequency (Hz)")
ylabel("$\psi$","Interpreter","Latex")

nexttile
plot(freq,real(W),'s',freq,imag(W),'^')
grid
legend("Real","Imag")
xlabel("Frequency (Hz)")
ylabel("W")
%Neglecting Willis Coupling------------------------------------------------
%Forward (R)
r_f = sqrt((1-R.*R+T.^2).^2 - T.^2);
x_f = (1-R.*R+T.^2+r)./(2.*T);

Z_f = (Z0.*r)./((1-R).*(1-R)-T.^2);
k_f = 1/L.*(1i.*log(x)+2*pi.*m);

kappa_f = Z_f.*omega./k_f;
rho_f = Z_f.*k_f./omega;

%Backward
r_b = sqrt((1-R_B.*R_B+T.^2).^2 - T.^2);
x_b = (1-R_B.*R_B+T.^2+r)./(2.*T);

Z_b = (Z0.*r)./((1-R_B).*(1-R_B)-T.^2);
k_b = 1/L.*(1i.*log(x)+2*pi.*m);

kappa_b = Z_b.*omega./k_b;
rho_b = Z_b.*k_b./omega;

f3 = figure();
figure(f3)
tiledlayout(1,2)
nexttile
plot(freq,real(rho_f),'s',freq,imag(rho_f),'^',freq,real(rho_b),'s',freq,imag(rho_b),'^')
grid
legend("Real Forward","Imag Forward","Real Backward","Imag Backward")
xlabel("Frequency (Hz)")
ylabel("$\rho$ ($kg/m^3$)","Interpreter","Latex")

nexttile
plot(freq,real(kappa_f)/1000,'s',freq,imag(kappa_f)/1000,'^',freq,real(kappa_b)/1000,'s',freq,imag(kappa_b)/1000,'^')
grid
legend("Real Forward","Imag Forward","Real Backward","Imag Backward")
xlabel("Frequency (Hz)")
ylabel("$\kappa$ (kPa)","Interpreter","Latex")