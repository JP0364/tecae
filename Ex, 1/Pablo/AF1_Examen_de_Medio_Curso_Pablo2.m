close all
clc
clear

%Pablo Arturo Balboa Sanchez 2077799
%Materiales Compuestos 2
%Examen Medio Curso
%%
%Problema 5

Ef1=230; %GPa
Ef2=25; %GPa
Em=3.3; %GPa
vm=0.35;
vf12=0.22;
vf21=vf12;

Vf=0.54; %Fraccion Volumetrica de la Fibra
Vm=1-Vf;

%Calculo de E1
E1 = (Ef1*Vf) + (Em*Vm);

%Calculo de E2
E2=(Vf/Ef2)+(Vm/Em);
E2=(E2)^-1;

%Calculo de Poisson 12
vc12=(vf12*Vf)+(vm*Vm);

fprintf('Resultados del Examen de Medio Curso\n');
fprintf('Pablo Arturo Balboa Sanchez\n');
fprintf('2077799\n');
fprintf('Resultados del Problema 5:\n');
% Imprimir el valor de E2 en la ventana de comandos
fprintf('E1: %.4f GPa\n', E1);
fprintf('E2: %.4f GPa\n', E2);
fprintf('Nu 12: %.4f GPa\n', vc12);

%%
%Problema 6
Gf12=28.3; %GPa
Gm=1.27; %GPa
Vf=0.57; %Vf correct to meet design specifications
Vm=1-Vf;

%ROM Equations
G12=(Vf/Gf12)+(Vm/Gm);
G12=(G12)^-1;
% Convertir E2 de Pascales a Gigapascales
fprintf('Resultados del problema 6:\n');
% Imprimir el valor de E2 en la ventana de comandos
fprintf('ROM G12: %.4f GPa\n', G12);

%Modify ROM Equation
n=0.6;
G12=((Vf/Gf12)+(n*Vm/Gm))/(Vf+(n*Vm));
G12=(G12)^-1;
fprintf('Modify ROM G12: %.4f GPa\n', G12);

%Elasticity Method
G12 = Gm * (( (Gm + Gf12) - Vf * (Gm - Gf12) ) / ( (Gm + Gf12) + Vf * (Gm - Gf12) ));
fprintf('Elasticity Method G12: %.4f GPa\n', G12);

%%
%Problema 7

%Sigma en el eje X en Pascales
SigmaX=230; %MPa
SigmaY=-50; %MPa
Tao12=-45; %MPa

% Módulos de Young para la matriz en Pascales
Young1 = 181000; %MPa
Young2 = 10300; %MPa
Young3=Young2;

% Ángulo Theta 
Theta = 0; % grados con incrementos de 1 grado

% Módulos de cortante para la matriz en Pascales
G12 = 7170; %MPa
G13 = G12;
%G23 = 3.28e9;

% Coeficientes de Poisson para la matriz
v12 = 0.28;
v13 = v12;
%v23 = 0.428;

L = 50e-3; % Lado de compuesto en metros

% Coeficientes de Poisson
v21 = (Young2 * v12) / Young1;
v31 = (Young3 * v13) / Young1;
%v32 = (Young3 * v23) / Young2;

% Matriz de cumplimiento S reducida
S_redox = [1/Young1, -v12/Young1, 0; 
           -v12/Young1, 1/Young2, 0;
           0, 0, 1/G12];

%Sigma vector
Sigma=[SigmaX SigmaY Tao12]';

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


fprintf('Respuestas para el problema 7:\n');
fprintf('Strain 1: %.4f\n', epsilon(1));
fprintf('Strain 2: %.4f\n', epsilon(2));
fprintf('Gamma12: %.4f\n', epsilon(3));

%%
%Problema 8
clear
%Arista
L = 50; %cm

%Deformación
Delta1=0.2329; %cm
Delta2=(-0.0291)*2; %cm

%Deformacion unitaria
epsilon1=Delta1/L;
epsilon2=Delta2/L;
Gamma12=0;

% Módulos de Young para la matriz en Pascales
Young1 = 70000; %MPa
Young2 = 4000; %MPa
Young3=Young2;

% Ángulo Theta 
Theta = 0; % grados con incrementos de 1 grado

% Módulos de cortante para la matriz en Pascales
G12 = 1;
G13 = G12;
%G23 = 3.28e9;

% Coeficientes de Poisson para la matriz
v12 = 0.25;
v13 = v12;
%v23 = 0.428;

% Coeficientes de Poisson
v21 = (Young2 * v12) / Young1;
v31 = (Young3 * v13) / Young1;
%v32 = (Young3 * v23) / Young2;

% Matriz de cumplimiento S reducida
S_redox = [1/Young1, -v12/Young1, 0; 
           -v12/Young1, 1/Young2, 0;
           0, 0, 1/G12];

%Matriz S inversa
C_redox=inv(S_redox);

%Vector de entrada de Strain
Strain=[epsilon1 epsilon2 Gamma12]';

%Vector Resultado
Stress=C_redox*Strain;

fprintf('Respuestas para el problema 8:\n');
fprintf('Stress 1: %.4f MPa\n', Stress(1));
fprintf('Stress 2: %.4f MPa\n', Stress(2));
fprintf('Tao12: %.4f MPa\n', Stress(3));

%%
%Problema 9

%Sigma en el eje X en Pascales
SigmaX=230; %MPa
SigmaY=0; %MPa
TaoXY=0; %MPa

% Módulos de Young para la matriz en MegaPascales
Young1 = 142*1000; %MPa
Young2 = 8.6*1000; %MPa
Young3=Young2;

% Ángulo Theta 
Theta = 20; % grados con incrementos de 1 grado

%Sigma vector
Sigma=[SigmaX SigmaY TaoXY]';

%Calculo de Matrix T
    m = cos(Theta * pi / 180);
    n = sin(Theta * pi / 180);
    T = [m*m, n*n, 2*m*n;
         n*n, m*m, -2*m*n;
         -m*n, m*n, m*m - n*n];

%Calculo de Esfuerzos en el sistema primario
Epsilon_primario=T*Sigma;

fprintf('Respuestas para el problema 9\n');
fprintf('Epsilon 1: %.4f MPa\n', Epsilon_primario(1));
fprintf('Epsilon 2: %.4f MPa\n', Epsilon_primario(2));
fprintf('Tao12: %.4f\n', Epsilon_primario(3));