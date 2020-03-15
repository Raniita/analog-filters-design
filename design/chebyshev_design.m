%% Filtro paso bajo Chebyshev

% Ganancia maxima banda pasante 20dB
% Frecuencia corte 1.2kHz
% Att min requerida -25dB a 1.7kHz
% Implementacion: etapas ganancia infinita y realimentacion multiple
% Rizado banda pasante 1 dB
% normalizar Wp,Ws. wc/w = 1/ws. Entonces, wc=1.2kHz, w=1.7kHz, 
% ws=1,416
Wp=1;
Ws=1.416;
Rp = 1;
Rs = 25;

% No conozco el orden del filtro, ergo cheb1ord
[N, Wp1] = cheb1ord(Wp, Ws, Rp, Rs,'s');

% Realmente Wp1 me va a coincidir con Wp, me interesa N
n=N;

% Polos usando Matlab
[z,p,k] = cheb1ap(n,Rp);
zplane(z,p);
title('Diagrama polos y ceros');

% Localizacion de polos
orden = 5;
Ap = 10^0.1;
epsilon = sqrt(10^(0.1*Ap)-1);
k = 0:1:orden-1;
sk = -sin((pi/2*orden)*(1+2*k))*sinh((1/orden)*asinh(1/epsilon)) + j*cos((pi/2*orden)*(1+2*k))*cosh((1/orden)*asinh(1/epsilon))

alpha = epsilon^-1 + sqrt((epsilon^(-2) )-1);
a = (alpha^(1/orden)-alpha^(-1/orden))/(2);
b = (alpha^(1/orden)+alpha^(-1/orden))/(2);

pk = -a*sin(((pi)/(2*orden))*(1+2*k)) + j*b*cos(((pi)/(2*orden))*(1+2*k))

% Agrupamos
c1 = abs(p(1) + p(5));
c2 = abs(p(2) + p(4));
c3 = abs(p(3));
c3 = [1 abs(p(3))];

% Ahora creo la función H(s) del filtro, tengo que tener en cuenta 
% que me va a devolver el filtro normalizado, wn no puede ser
% mayor a 1 en esta funcion
wn = 1;
n = N;

[b,a] = cheby1(n,Rp,wn,'s');

p1 = tf(b,a);

% Other way to go to the same place 
% [A,B,C,D] = zp2ss(z,p,k);
% [At,Bt,Ct,Dt] = lp2lp(A,B,C,D,1200*pi*2);
% [bc,ac] = ss2tf(At,Bt,Ct,Dt);       % Convert to TF form.

% Ahora vamos a ajustarlo a la frec de corte que queremos
% wc tiene que estar en rad/sec

wc = 1200*2*pi;

[bt,at] = lp2lp(b,a,wc);

p2=tf(bt,at);

% Representacion en frecuencia

Wpoints = 0.01:0.025:10000; %Para una resolucion de 0.01rad/s

%h = bodeplot(p1,'b',p2,'r',Wpoints);
%P = getoptions(h);
%P.FreqUnits = 'rad/s';
%P.MagUnits = 'db';
%P.Title.String = 'Chebyshev filters';
%P.PhaseVisible = 'off';
%P.XLimMode = 'manual';
%P.XLim = ([0.01 10]);
%P.YLimMode = 'manual';
%P.YLim = ([-4 1]);
%P.Grid = 'on';
%setoptions(h,P);

P = bodeoptions;
P.FreqUnits = 'rad/s';
P.MagUnits = 'db';
P.Title.String = 'Filtro Chebyshev puesto B5';
P.PhaseVisible = 'off';
P.XLimMode = 'manual';
P.XLim = ([0.01 15000]);
P.YLimMode = 'manual';
P.YLim = ([-4 1]);
P.Grid = 'on';
h = bodeplot(p1,'b',p2,'r',Wpoints,P)

% Descomponemos en polos y residuos
[residuos, polos, k] = residue(bt,at);

% Calculamos los numeradores
num1 = abs((residuos(3)^2)*(residuos(4)^2));
num2 = abs((residuos(1)^2)*(residuos(2)^2));
num3 = residuos(5);

% Calculamos los denominadores
den1 = [1 -(polos(3) + polos(4)) (polos(3)^2)*(polos(4)^2)];
den2 = [1 -(polos(1) + polos(2)) (polos(1)^2)*(polos(2)^2)];
den3 = [1 polos(5)];

% Primer cuadratico
wo1 = sqrt(den1(3));
alpha1 = den1(2)/wo1;
Q1 = alpha1^(-1);
Ho1 = num1/den1(3);

% Segundo cuadratico
wo2 = sqrt(den2(3));
alpha2 = den2(2)/wo2;
Q2 = alpha2^(-1);
Ho2 = num2/den2(3);

% Simple
wo3 = den3(2);
Q3 = 1;
Ho3 = num3/wo3;

% Etapas con >Q van antes.
%Ho = at(6)
K = (10^(20/20)) * at(6);
K_norm = 1 * at(6);
K_calc = num1 + num2 + num3;

% Representamos las etapas
p_cuadratico1 = tf(num1, den1);
p_cuadratico2 = tf(num2, den2);
p_simple = tf(num3, den3);

bode(p_cuadratico1)
bode(p_cuadratico2)
bode(p_simple)

% Utilizar tiff print(tiff) 300px
%print -dtiff -r300 ./test2