close all
clc
clear

%Materiales Compuestos II
%AF 6 Programación de teorias de laminados
%Pablo Arturo Balboa Sanchez Mat. 2077799
%%Problema 1
%Inciso a)

%Sigma en el eje X en Pascales
SigmaX=100e6;

% Módulos de Young para la matriz en Pascales
Young1 = 50.0e9;
Young2 = 15.20e9;
Young3=Young2;

% Ángulo Theta 
Theta = 0; % grados con incrementos de 1 grado

% Módulos de cortante para la matriz en Pascales
G12 = 4.70e9;
G13 = G12;
G23 = 3.28e9;

% Coeficientes de Poisson para la matriz
v12 = 0.254;
v13 = v12;
v23 = 0.428;

L = 50e-3; % Lado de compuesto en metros

% Coeficientes de Poisson
v21 = (Young2 * v12) / Young1;
v31 = (Young3 * v13) / Young1;
v32 = (Young3 * v23) / Young2;

% Matriz de cumplimiento S reducida
S_redox = [1/Young1, -v12/Young1, 0; 
           -v12/Young1, 1/Young2, 0;
           0, 0, 1/G12];

%Sigma vector
Sigma=[SigmaX 0 0]';

%Calculo de Sbar
    m = cos(Theta * pi / 180);
    n = sin(Theta * pi / 180);
    T = [m*m, n*n, 2*m*n;
         n*n, m*m, -2*m*n;
         -m*n, m*n, m*m - n*n];
    Tinv = inv(T);
    Sbar = Tinv * S_redox * T;
    
%Deformacion Unitaria
epsilon=Sbar*Sigma;

%Deformacion unitaria
deltax=(L)*epsilon(1);
deltay=(L)*epsilon(2);


%Deformacion total en metros
dx=deltax+L;
dy=deltay+L;
gammaxy=epsilon(3);

%Distancia en X deformada en milimetros
dx=dx/1e-3;
dy=dy/1e-3;

fprintf('Respuestas para el problema # 1 de la AF 6\n');
fprintf('a) Dimensiones deformadas para Theta (\x03B8) = 0°\n');
fprintf('DeltaX (\x0394X): %.4f mm\n', dx);
fprintf('DeltaY (\x0394Y): %.4f mm\n', dy);
fprintf('GammaXY (\x03B3XY): %.4f\n', gammaxy);

%Inciso b)
%Reestableciendo Theta
Theta=45;

%Calculo de Sbar
    m = cos(Theta * pi / 180);
    n = sin(Theta * pi / 180);
    T = [m*m, n*n, 2*m*n;
         n*n, m*m, -2*m*n;
         -m*n, m*n, m*m - n*n];
    Tinv = inv(T);
    Sbar = Tinv * S_redox * T;
    
%Deformacion Unitaria
epsilon=Sbar*Sigma;

%Deformacion unitaria
deltax=(L)*epsilon(1);
deltay=(L)*epsilon(2);

%Deformacion total en metros
dx=deltax+L;
dy=deltay+L;
gammaxy=epsilon(3);

%Distancia en X deformada en milimetros
dx=dx/1e-3;
dy=dy/1e-3;

fprintf('b) Dimensiones deformadas para Theta (\x03B8) = +45°\n');
fprintf('DeltaX (\x0394X): %.4f mm\n', dx);
fprintf('DeltaY (\x0394Y): %.4f mm\n', dy);
fprintf('GammaXY (\x03B3XY): %.4f\n', gammaxy);

%Inciso c)
%Reestableciendo Theta
Theta=-45;

%Calculo de Sbar
    m = cos(Theta * pi / 180);
    n = sin(Theta * pi / 180);
    T = [m*m, n*n, 2*m*n;
         n*n, m*m, -2*m*n;
         -m*n, m*n, m*m - n*n];
    Tinv = inv(T);
    Sbar = Tinv * S_redox * T;
    
%Deformacion Unitaria
epsilon=Sbar*Sigma;

%Deformacion en metros
deltax=(L)*epsilon(1);
deltay=(L)*epsilon(2);


%Deformacion total en metros
dx=deltax+L;
dy=deltay+L;
gammaxy=epsilon(3);

%Distancia en X deformada en milimetros
dx=dx/1e-3;
dy=dy/1e-3;

fprintf('c) Dimensiones deformadas para Theta (\x03B8) = -45°\n');
fprintf('DeltaX (\x0394X): %.4f mm\n', dx);
fprintf('DeltaY (\x0394Y): %.4f mm\n', dy);
fprintf('GammaXY (\x03B3XY): %.4f\n', gammaxy);

%%Segundo problema de la AF 6

% Coeficientes de Poisson
%v21 = (Young2 * v12) / Young1;
%v31 = (Young3 * v13) / Young1;
%v32 = (Young3 * v23) / Young2;
v21 = v12; % Asumiendo que el material se comporte isotropicamente 

% Rango de Theta
theta_grados = -90:10:90;  % grados

% Creacion de matrices para los valores de las constantes
Youngx = zeros(size(theta_grados)); %Modulo Ex
Youngy = zeros(size(theta_grados)); %Modulo Ey
Gxy = zeros(size(theta_grados)); %Modulo de Corte xy
Poisson_xy = zeros(size(theta_grados)); % Poisson xy
Poisson_yx = zeros(size(theta_grados)); % Poisson yx

% Calculo de las constantes 
for i = 1:length(theta_grados)
    % Indice para actualizar el valor de Theta en cada ciclo
    theta = theta_grados(i);
    
    % Calculo de m y n, pasando de grados a radianes
    m = cos(theta* pi / 180);
    n = sin(theta* pi / 180);
    
    % Calculo para el modulo de Poisson xy
    den_Poisson_xy = m^4 + (Young1 / G12 - 2 * v12) * n^2 * m^2 + (Young1 / Young2) * n^2;
    num_Pisson_xy = v12 * (n^4 + m^4) - (1 + Young1 / Young2 - Young1 / G12) * n^2 * m^2;
    Poisson_xy(i) = num_Pisson_xy / den_Poisson_xy;
    
    % Calculo para el modulo de Young Ey
    den_Ey = m^4 + (Young2 / G12 - 2 * v21) * n^2 * m^2 + (Young2 / Young1) * n^4;
    Youngy(i) = Young2 / den_Ey;
    
    % Calculo para modulo de corte Gxy
    den_Gxy = n^4 + m^4 + 2 * (2 * G12 * (1 + 2 * v12) / Young1 + 2 * G12 / Young2 - 1) * n^2 * m^2;
    Gxy(i) = G12 / den_Gxy;
    
    % Calculo de modulo de Young Ex
    denom_Ex = m^4 + (Young1 / G12 - 2 * v12) * n^2 * m^2 + (Young1 / Young2) * n^4;
    Youngx(i) = Young1 / denom_Ex;
    
    % Calculo del modulo de Poisson yx
    denom_Vyx = m^4 + (Young2 / G12 - 2 * v21) * n^2 * m^2 + (Young2 / Young1) * n^2;
    num_Vyx = v21 * (n^4 + m^4) - (1 + Young2 / Young1 - Young2 / G12) * n^2 * m^2;
    Poisson_yx(i) = num_Vyx / denom_Vyx;
end

% Convertir las constantes a GPa para la gráfica
Youngx = Youngx / 1e9;
Youngy = Youngy / 1e9;
Gxy = Gxy / 1e9;

% Graficar las constantes elásticas en figuras separadas

% Gráfica de Ex
figure;
plot(theta_grados, Youngx, 'r', 'LineWidth', 2);
grid on;
xlabel('\theta (grados)');
ylabel('E_x (GPa)');
title('Módulo de Young E_x en función de \theta');

% Gráfica de Ey
figure;
plot(theta_grados, Youngy, 'g', 'LineWidth', 2);
grid on;
xlabel('\theta (grados)');
ylabel('E_y (GPa)');
title('Módulo de Young E_y en función de \theta');

% Gráfica de Gxy
figure;
plot(theta_grados, Gxy, 'b', 'LineWidth', 2);
grid on;
xlabel('\theta (grados)');
ylabel('G_{xy} (GPa)');
title('Módulo de Corte G_{xy} en función de \theta');

% Gráfica de Poisson xy
figure;
plot(theta_grados, Poisson_xy, 'm', 'LineWidth', 2);
grid on;
xlabel('\theta (grados)');
ylabel('\nu_{xy}');
title('Módulo de Poisson \nu_{xy} en función de \theta');

% Gráfica de Poisson yx 
figure;
plot(theta_grados, Poisson_yx, 'k', 'LineWidth', 2);
grid on;
xlabel('\theta (grados)');
ylabel('\nu_{yx}');
title('Módulo de Poisson \nu_{yx} en función de \theta');