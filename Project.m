%clear all;
%close all;
format short;

syms l mu E r_o s A I nus nu;
syms x w w1 w2 w3 P P1 P2 P3;
syms Jf Jt mf mt M L;
syms n a b K_a J JBar L_1 L_2 A1 A2 beta ep v_o;

syms omega1 omega2;

s = 10;
l = 0.1*(10^(-9));
nu = 0.23;
v_o = 1000;

n = 5;
a = 0.22;
r_o = 0.87*(10^(-9));
E = 10^(12);
A = 1.495398*(10^(-18));
nus = (1-nu)/((1+nu)*(1-2*nu));  %1.097858;
mu = 4.2*(10^(11));
I = 3.879810*(10^(-37));
mf = 8.824734*(10^(-16));
mt = 3.439415*(10^(-15));
M = mf + mt;
L = 2*r_o*s;

Jf = (1/12)*mf*(10^(-18))*( 0.8427 + 3.0276*s*s )
Jt = (1/12)*mt*(10^(-18))*( 3.1134 + 3.0276*s*s )


b = (I*nus)/(4*A*r_o*r_o*s*s)
K_a = (l*l*mu)/(16*E*r_o*r_o*s*s)
J = (Jf + Jt);
JBar = J/(M*L*L)

L_1 = abs(sqrt(((((2*n+1)*pi)^2)*(K_a+b)/4 - a)/JBar))
L_2 = L_1;
beta = abs(power((L_2*L_2)/K_a, 0.25))

A1 = -(cosh(beta)+cos(beta))/(sqrt(2)*beta*sinh(beta)*sin(beta));
A2 = (sinh(beta)+sin(beta))/(sqrt(2)*beta*sinh(beta)*sin(beta));
%A1 = roundn(A1, -6)
%A2 = roundn(A2, -6)

w = A2*(cosh(beta*x) - cos(beta*x)) + A1*(sinh(beta*x) - sin(beta*x));
w1 = diff(w);
w2 = diff(w, 2);
w3 = diff(w, 3);
w4 = diff(w, 4);

P = sqrt(2)*sin((2*n+1)*pi*x/2);
P1 = diff(P);
P2 = diff(P, 2);
P3 = diff(P, 3);

syms C1 C2 C3 C4 C5 C6 C7 CS Kappa;
 
C1 = eval(int(w*w, 0, 1))
C2 = a*eval(int(w*w2, 0, 1)) + K_a*eval(int(w*w4, 0, 1))
C3 = eval(int(w*w2, 0, 1))
C4 = a*eval(int(w*P1 + K_a*w*P3, 0, 1)) - K_a*eval(int(w*P3, 0, 1))
C5 = JBar*eval(int(P*P, 0, 1))
C6 = a*eval(int(P*w1 - K_a*P*w3, 0, 1)) - K_a*eval(int(P*w3, 0, 1))
C7 = (b+K_a)*eval(int(P*P2 + a*P*P, 0, 1)) + a*eval(int(P*P, 0, 1))

Kappa = abs(2*sqrt(mf/M))
CS = Kappa*eval(int(w*w1, 0, 1))

syms a1 a2 a3 a4 a1s a2s a2C a4C;
a1 = -C7/C5 - C2/C1
a2 = C3/C1
a3 = (C2*C7 - C4*C6)/(C1*C5)
a4 = -(C3*C7)/(C1*C5)
a1s = CS/C1
a2s = -(CS*C7)/(C1*C5)
ep = 0.1;
a2C = a2/ep;
a4C = a4/ep;

omega1 = abs(sqrt( a1/2 - sqrt(a1*a1/4 - a3) ))
omega2 = abs(sqrt( a1/2 + sqrt(a1*a1/4 - a3) ))

syms Rf1 Rf2 Rf3 Rf4 Rf5 Rf6 Rf7 Rf8
Rf1 = (omega2 - omega1)/2;
Rf2 = (omega2 - omega1);
Rf3 = (omega2 + omega1)/2;
Rf4 = (omega2 + omega1);
Rf5 = omega1;
Rf6 = 2*omega1;
Rf7 = omega2;
Rf8 = 2*omega2;

syms sigma OmegaBar
syms b1 b2 b3 v_D
syms Tao Chi_1 Chi_2P Y1 Y2

Tao = (omega1*omega1 - a4C/a2C)/(omega1*(2*omega1*omega1 - a1));

syms X1 X2 X3 X4 X5 X6
  
OmegaBar = 2.*omega1;
X1 = (omega1-2*OmegaBar)^2/((omega1-2*OmegaBar)^4 - a1*(omega1-2*OmegaBar)^2 + a3);
X2 = (omega1+2*OmegaBar)^2/((omega1+2*OmegaBar)^4 - a1*(omega1+2*OmegaBar)^2 + a3);
X3 = (a2C*(omega1+OmegaBar)^2 - a4C)/((omega1+OmegaBar)^4 - a1*(omega1+OmegaBar)^2 + a3);
X4 = 1/((omega1+2*OmegaBar)^4 - a1*(omega1+2*OmegaBar)^2 + a3);
X5 = 1/((omega1-2*OmegaBar)^4 - a1*(omega1-2*OmegaBar)^2 + a3);
X6 = (a2C*(omega1-2*OmegaBar)^2 - a4C)/((omega1-2*OmegaBar)^4 - a1*(omega1-2*OmegaBar)^2 + a3);

XA = [];
XB = [];
Y = [];

for v_D = 0:0.01:500

    b1 = a2C*(mf/(E*A))*(v_o^2 + v_D^2/2);
    b2 = (a2C/2)*(mf/(E*A))*(v_D^2);
    b3 = (2*a2C)*(mf/(E*A))*(v_o*v_D);
    
    Chi_1 = a2C*Tao*(b2*b2/2)*(X1 + X2) + Tao*(b3*b3/4)*X3 - a4C*Tao*(b2*b2/4)*(X4 + X5) + (omega1*a2C*(b1*b1 + b3*b3/4) + (6*omega1^2 - a1)*(b3*b3/4-b1*b1)*Tao/4)*Tao/(omega1*(2*omega1^2-a1));
    Chi_2P = Tao*(b2*b3/4)*(X3 + X6) + omega1*a2C*b1*b3*Tao;
    
    Y1 = (b1*Tao + ep*Chi_1 - b3*Tao/2 - ep*Chi_2P)/(2 - ep*(6*omega1*omega1-a1)*b3*Tao/2);
    Y2 = (b1*Tao + ep*Chi_1 + b3*Tao/2 + ep*Chi_2P)/(2 + ep*(6*omega1*omega1-a1)*b3*Tao/2);
    
    XA = [XA, 2*omega1+2*ep*Y1];
    XB = [XB, 2*omega1+2*ep*Y2];
    Y = [Y, v_D];
end 


xlabel('\Omega')
ylabel('V_D');

plot(XA, Y, ':k');
hold on
plot(XB, Y, ':k');
hold on