%% Filtro paso bajo Chebyshev
% Ganancia maxima banda pasante 20dB
% Frecuencia corte 1.2kHz
% Att min requerida -25dB a 1.7kHz
% Implementacion: etapas ganancia infinita y realimentacion multiple
% Rizado banda pasante 1 dB

%% Calculo parametros primarios: epsilon, att min requerida, orden

r_db = 1;
r_lineal = 10^(r_db/10);
epsilon = sqrt(r_lineal-1);
gamma = -25; %atenuacion minima requerida

A = 10^(abs(gamma)/20);

w = 1.7*10^3; %valor al que se da A
wc = 1.2*10^3; %valor corte filtro

ws = w/wc;

g = sqrt((A^2-1)/(epsilon^2)); %parametro que me ayuda a calcular nmin

nmin = log(g+sqrt(g^2-1))/log(ws+sqrt(ws^2-1));

n = round(nmin); %el orden del filtro no puede ser decimal

% Comprobacion matlab

Wp=1;
Ws=1.416;
Rp = 1;
Rs = 25;

[N, Wp1] = cheb1ord(Wp, Ws, Rp, Rs,'s');

%% Calculo polos y denominadores
alpha = epsilon^(-1) + sqrt(epsilon^(-2)+1);

a = (alpha^(1/n)- alpha^(-1/n))/2;

b = (alpha^(1/n) + alpha^(-1/n))/2;

k_pk = 0:1:n-1;

pk = -a*sin((pi/10)*(1+2*k_pk)) + i*b*cos((pi/10)*(1+2*k_pk));

c1 = [1 abs(sum(pk(1)+pk(5))) abs(imag(pk(1))*imag(pk(5)))];
c2 = [1 abs(sum(pk(2)+pk(4))) abs(imag(pk(2))*imag(pk(4)))];
l1 = [1 abs(pk(3))];

% Comprobacion matlab
n=N;

[z,p,k] = cheb1ap(n,Rp);

zplane(z,p)
%% Construccion funciones transferencia

%Normalizado
b_lineal = [0.2886];
a_lineal = [1 0.2886];

b_cuadratico1 = [0.9883];
a_cuadratico1 = [1 0.1789 0.9883];

b_cuadratico2 = [0.4293];
a_cuadratico2 = [1 0.4684 0.4293];

% Desplazado en freq a 1.2kHz
bt_lineal = [2175.99];
at_lineal = [1 2175.99];

bt_cuadratico1 = [56*10^6];
at_cuadratico1 = [1 1348.87 56*10^6];

bt_cuadratico2 = [24.4*10^6];
at_cuadratico2 = [1 3531.7 24.4*10^6];

% Obtenemos las funciones de transferencia
p_lineal = tf(b_lineal, a_lineal);
p_cuadra1 = tf(b_cuadratico1, a_cuadratico1);
p_cuadra2 = tf(b_cuadratico2, a_cuadratico2);

pt_lineal = tf(bt_lineal, at_lineal);
pt_cuadra1 = tf(bt_cuadratico1, at_cuadratico1);
pt_cuadra2 = tf(bt_cuadratico2, at_cuadratico2);

pt_total = pt_lineal * pt_cuadra1 * pt_cuadra2;

% Comprobacion matlab

% tengo que tener en cuenta que me va a devolver el filtro normalizado, 
% wn no puede ser mayor a 1 en esta funcion

wn = 1;
n = N;

[b,a] = cheby1(n,Rp,wn,'s');

p1 = tf(b,a);

% Ahora vamos a ajustarlo a la frec de corte que queremos
% wc tiene que estar en rad/sec


wc = 1200*2*pi;

[bt,at] = lp2lp(b,a,wc);

p2=tf(bt,at);

%% Representacion en frecuencia
%Wpoints = 0.01:0.05:20000; %Para una resolucion de 0.01rad/s

P = bodeoptions;
P.FreqUnits = 'rad/s';
P.MagUnits = 'db';
P.Title.String = 'Filtro Chebyshev puesto B5';
%P.PhaseVisible = 'off';
P.XLimMode = 'auto';
%P.XLim = ([0.01 15000]);
P.YLimMode = 'auto';
%P.YLim = ([-4 20]);
P.Grid = 'on';

%bodeplot(pt_total, 'b', P);
bodeplot(pt_lineal, 'b', pt_cuadra1, 'r', pt_cuadra2, 'g', P)
%bodeplot(p_lineal, 'b', p_cuadra1, 'r', p_cuadra2, 'g', pt_lineal, 'b', pt_cuadra1, 'r', pt_cuadra2, 'g', P)

%% Sintesis electronica
% Usando ganancia infinita y realimentacion multiple
% Ordenamos segun Q> primero

% Etapa cuadratica 1. Q=5.55
Ho1 = bt_cuadratico1/at_cuadratico1(3);
wo1 = sqrt(at_cuadratico1(3));
Q1 = wo1/at_cuadratico1(2);
alpha1 = Q1^(-1);
refC1 = 1*10^-9;        % 1nF

[R1_cuadr1, C2_cuadr1, R3_cuadr1, R4_cuadr1, C5_cuadr1] = lp_reaMult(Ho1, alpha1, wo1, refC1);

% Etapa cuadratica 2. Q=1.39
Ho2 = bt_cuadratico2/at_cuadratico2(3);
wo2 = sqrt(at_cuadratico2(3));
Q2 = wo2/at_cuadratico2(2);
alpha2 = Q2^(-1);
refC2 = 2.2*10^-9;     % 2.2nF

[R1_cuadr2, C2_cuadr2, R3_cuadr2, R4_cuadr2, C5_cuadr2] = lp_reaMult(Ho2, alpha2, wo2, refC2);

% Etapa Lineal. Filtro RC
wo3 = 2175.99;
refC3 = 470*10^-9;      % 470nF

R_lineal = 1/(wo3*refC3);
