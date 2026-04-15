clear; clc; close all;

%% =========================================================
% ITEM 1 - MODELO RLC EN ESPACIO DE ESTADOS
%% =========================================================

% Parametros del circuito
R = 2200;        % ohm
L = 500e-3;      % H
C = 10e-6;       % F

% Matrices de estado
A = [-R/L   -1/L;
      1/C    0];

B = [1/L;
     0];

% Salida: caida en la resistencia
C_mat = [R 0];

D = 0;

% Sistema
sys = ss(A,B,C_mat,D);

% Tiempo
t = 0:1e-4:2;

% Entrada escalon alternante
u = 12 * sign(sin(2*pi*t));

% Simulacion
[y,t,x] = lsim(sys,u,t);

% Graficas ITEM 1
figure

subplot(4,1,1)
plot(t,x(:,1),'LineWidth',1.5)
grid on
ylabel('i(t) [A]')
title('Corriente')

subplot(4,1,2)
plot(t,x(:,2),'LineWidth',1.5)
grid on
ylabel('V_c(t) [V]')
title('Tension en el capacitor')

subplot(4,1,3)
plot(t,u,'LineWidth',1.5)
grid on
ylabel('V_e(t) [V]')
title('Entrada')

subplot(4,1,4)
plot(t,y,'LineWidth',1.5)
grid on
ylabel('V_r(t) [V]')
xlabel('Tiempo [s]')
title('Salida')


%% =========================================================
% ITEM 2 - DATOS EXPERIMENTALES + METODO DE CHEN
%% =========================================================

% Cargar datos del Excel
data = readmatrix('Curvas_Medidas_RLC_2026.xls');

t_data  = data(:,1);
i_data  = data(:,2);
vc_data = data(:,3);
ve_data = data(:,4);

% Graficas de los datos
figure

subplot(3,1,1)
plot(t_data,ve_data,'LineWidth',1.2)
grid on
ylabel('V_e(t)')
title('Entrada medida')

subplot(3,1,2)
plot(t_data,vc_data,'LineWidth',1.2)
grid on
ylabel('V_c(t)')
title('Tension en el capacitor')

subplot(3,1,3)
plot(t_data,i_data,'LineWidth',1.2)
grid on
ylabel('i(t)')
xlabel('Tiempo [s]')
title('Corriente')


%% =========================================================
% SELECCION DEL ESCALON
%% =========================================================

idx = t_data >= 0.1 & t_data <= 0.5;

t_step  = t_data(idx);
vc_step = vc_data(idx);
i_step  = i_data(idx);

% Tiempo relativo
t0 = t_step - t_step(1);

% Normalizacion
K = max(vc_step);
yn = vc_step / K;

figure
plot(t0,yn,'LineWidth',1.5)
grid on
title('Salida normalizada')
xlabel('Tiempo relativo')
ylabel('y_n(t)')


%% =========================================================
% METODO DE CHEN (SIN CERO)
%% =========================================================

% Elegimos t1
t1 = 0.05;

t_ch = [t1 2*t1 3*t1];

% Interpolacion
y_ch = interp1(t0,yn,t_ch);

y1 = y_ch(1);
y2 = y_ch(2);
y3 = y_ch(3);

fprintf('\nPuntos de Chen:\n')
fprintf('y(t1)  = %.6f\n',y1)
fprintf('y(2t1) = %.6f\n',y2)
fprintf('y(3t1) = %.6f\n',y3)

% Parametros k
k1 = y1 - 1;
k2 = y2 - 1;
k3 = y3 - 1;

% Calculo de b
b = 4*k1^3*k3 - 3*k1^2*k2^2 - 4*k2^3 + k3^2 + 6*k1*k2*k3;

% Calculo de alpha
alpha1 = (k1*k2 + k3 - sqrt(b)) / (2*(k1^2 + k2));
alpha2 = (k1*k2 + k3 + sqrt(b)) / (2*(k1^2 + k2));

% Constantes de tiempo
T1 = -t1 / log(alpha1);
T2 = -t1 / log(alpha2);

fprintf('\nResultados Chen:\n')
fprintf('T1 = %.6f s\n',T1)
fprintf('T2 = %.6f s\n',T2)

%% =========================================================
% MODELO IDENTIFICADO
%% =========================================================

s = tf('s');

G_chen = 1 / ((T1*s + 1)*(T2*s + 1));

% Respuesta del modelo
y_model = step(G_chen,t0);

figure
plot(t0,yn,'b','LineWidth',1.8); hold on
plot(t0,y_model,'k--','LineWidth',1.8)
plot(t_ch,y_ch,'ro','MarkerSize',8,'LineWidth',1.5)
grid on
legend('Medida','Chen','Puntos Chen')
title('Comparacion modelo Chen')
xlabel('Tiempo')
ylabel('Salida')


%% =========================================================
% OBTENCION DE R, L y C
%% =========================================================

% Coeficientes
LC = T1 * T2;
RC = T1 + T2;

fprintf('\nCoeficientes:\n')
fprintf('LC = %.8f\n',LC)
fprintf('RC = %.6f\n',RC)

% Derivada numerica
dvc = gradient(vc_step,t_step);

% Calculo de C
C_est = i_step ./ dvc;

% Filtrado
idx_valid = abs(dvc) > 0.1 * max(abs(dvc));

C_med = mean(C_est(idx_valid));

% Calculo de R y L
R_est = RC / C_med;
L_est = LC / C_med;

fprintf('\nParametros estimados:\n')
fprintf('C = %.6e F\n',C_med)
fprintf('R = %.3f Ohm\n',R_est)
fprintf('L = %.6f H\n',L_est)
