%% Design filtro paso bajo Chebyshev
% Ganancia maxima banda pasante 20dB
% Frecuencia corte 1.2kHz
% Att min requerida -25dB a 1.7kHz
% Implementacion: etapas ganancia infinita y realimentacion multiple
% Rizado banda pasante 1 dB

% Normalizado
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

% Representacion en frecuencia
%Wpoints = 0.01:0.05:20000; %Para una resolucion de 0.01rad/s

P = bodeoptions;
P.FreqUnits = 'rad/s';
P.MagUnits = 'db';
P.MagScale = 'linear';
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
% Usando ganancia infinita y realimentación múltiple
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
