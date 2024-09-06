% close all
clc
clear

%Pablo Arturo Balboa Sanchez 2077799
%Materiales Compuestos 2
%Actividad Fundamental 5
%%
%Problema 1

Ef1=85.6*10^9;
Ef2=14.8*10^9;
Em=3.45*10^9;
vm=0.36;
vf12=0.3;
vf21=vf12;
n=0.5;

Vf=0.65;
Vm=1-Vf;

%Inciso A)
E2=(Vf/Ef2)+(Vm/Em);
E2=(E2)^-1;
% Convertir E2 de Pascales a Gigapascales
E2_GPa = E2 / 10^9;

fprintf('Resultados de la Actividad Fundamental 5 Problema 1\n');
% Imprimir el valor de E2 en la ventana de comandos
fprintf('A) E2: %.4f GPa\n', E2_GPa);

%Inciso B)
E2=((Vf/Ef2)+(n*Vm/Em))/(Vf+(n*Vm));
E2=(E2)^-1;
E2_GPa = E2 / 10^9;

% Imprimir el valor de E2 en la ventana de comandos
fprintf('B) E2: %.4f GPa\n', E2_GPa);

%Inciso C)
nf = (Ef1*Vf + ((1 - vf12*vf21)*Em + vm*vf21*Ef1)*Vm) / (Ef1*Vf + Em*Vm);
nm = (((1 - vm^2)*Ef1 - (1 - vm*vf12)*Em)*Vf + Em*Vm) / (Ef1*Vf + Em*Vm);
E2=((nf*Vf)/Ef2)+((nm*Vm)/Em);
E2=(E2)^-1;
E2_GPa = E2 / 10^9;
fprintf('C) E2: %.4f GPa\n', E2_GPa);

%%
% Problema 2
%Inciso A
% Rango de Vf
Vf_range = 0:0.0005:1;

% Vector para datos de E2
E1_A_GPa = zeros(size(Vf_range));
E2_A_GPa = zeros(size(Vf_range));
E2_B_GPa = zeros(size(Vf_range));
E2_C_GPa = zeros(size(Vf_range));
E2_D_GPa = zeros(size(Vf_range));
E2_E_GPa = zeros(size(Vf_range));

% Ciclo para calcular E2 para cada Inciso
for i = 1:length(Vf_range)
    Vf = Vf_range(i);
    Vm = 1 - Vf;

    %Calculo de E1
    E1 = (Ef1*Vf) + (Em*Vm);
    E1_A_GPa(i) = E1 / 10^9;

    % Inciso A: Regla simple de mezclas
    E2_A = (Vf/Ef2) + (Vm/Em);
    E2_A = (E2_A)^-1;
    E2_A_GPa(i) = E2_A / 10^9;
    
    % Inciso B: Regla modificada con n = 0.4
    n = 0.4;
    E2_B = ((Vf/Ef2) + (n*Vm/Em)) / (Vf + (n*Vm));
    E2_B = (E2_B)^-1;
    E2_B_GPa(i) = E2_B / 10^9;
    
    % Inciso C: Regla modificada con n = 0.5
    n = 0.5;
    E2_C = ((Vf/Ef2) + (n*Vm/Em)) / (Vf + (n*Vm));
    E2_C = (E2_C)^-1;
    E2_C_GPa(i) = E2_C / 10^9;
    
    % Inciso D: Regla modificada con n = 0.6
    n = 0.6;
    E2_D = ((Vf/Ef2) + (n*Vm/Em)) / (Vf + (n*Vm));
    E2_D = (E2_D)^-1;
    E2_D_GPa(i) = E2_D / 10^9;
    
    % Inciso E: Usando nf y nm
    nf = (Ef1*Vf + ((1 - vf12*vf21)*Em + vm*vf21*Ef1)*Vm) / (Ef1*Vf + Em*Vm);
    nm = (((1 - vm^2)*Ef1 - (1 - vm*vf12)*Em)*Vf + Em*Vm) / (Ef1*Vf + Em*Vm);
    E2_E = ((nf*Vf)/Ef2) + ((nm*Vm)/Em);
    E2_E = (E2_E)^-1;
    E2_E_GPa(i) = E2_E / 10^9;
end

% Graficar todos los resultados
figure('Name', 'Problema_2');
plot(Vf_range, E2_A_GPa, '-o', 'LineWidth', 2, 'MarkerSize', 4,'DisplayName','Inciso A');
hold on;
plot(Vf_range, E2_B_GPa, '--', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Inciso B');
plot(Vf_range, E2_C_GPa, '-*', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Inciso C');
plot(Vf_range, E2_D_GPa, '*-', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Inciso D');
plot(Vf_range, E2_E_GPa, '-o', 'LineWidth', 2, 'MarkerSize', 6,'DisplayName','Inciso E');
% plot(Vf_range, E1_A_GPa, '-o', 'LineWidth', 2, 'MarkerSize', 6,'DisplayName','E1');
xlabel('Vf (Fracción Volumétrica de la Fibra)');
ylabel('E2 (GPa)');
title('Gráfico de E2 vs Vf (Problema 2)');
grid on;
legend;
hold off;

%%
%Problema 3
Gf12=28.3*10^9;
Gm=1.27*10^9;
Vf=0.55;
Vm=1-Vf;

%Inciso A
G12=(Vf/Gf12)+(Vm/Gm);
G12=(G12)^-1;
% Convertir E2 de Pascales a Gigapascales
G12_GPa = G12 / 10^9;
fprintf('Resultados de la Actividad Fundamental 5 Problema 3\n');
% Imprimir el valor de E2 en la ventana de comandos
fprintf('A) G12: %.4f GPa\n', G12_GPa);

%Inciso B
n=0.6;
G12=((Vf/Gf12)+(n*Vm/Gm))/(Vf+(n*Vm));
G12=(G12)^-1;
G12_GPa = G12 / 10^9;
fprintf('B) G12: %.4f GPa\n', G12_GPa);

%Inciso C
G12 = Gm * (( (Gm + Gf12) - Vf * (Gm - Gf12) ) / ( (Gm + Gf12) + Vf * (Gm - Gf12) ));
G12_GPa = G12 / 10^9;
fprintf('C) G12: %.4f GPa\n', G12_GPa);

%%
%Problema 4
Gf12=28.3*10^9;
Gm=1.27*10^9;
Vf_range = 0:0.1:1;

% Vector para datos de G12
G12_A_GPa = zeros(size(Vf_range));
G12_B_GPa = zeros(size(Vf_range));
G12_C_GPa = zeros(size(Vf_range));

% Ciclo para calcular E2 para cada Inciso
for i = 1:length(Vf_range)
    Vf = Vf_range(i);
    Vm = 1 - Vf;

%Inciso A
G12=(Vf/Gf12)+(Vm/Gm);
G12=(G12)^-1;
% Convertir E2 de Pascales a Gigapascales
G12_A_GPa(i) = G12 / 10^9;

%Inciso B
n=0.6;
G12=((Vf/Gf12)+(n*Vm/Gm))/(Vf+(n*Vm));
G12=(G12)^-1;
G12_B_GPa(i) = G12 / 10^9;

%Inciso C
G12 = Gm * (( (Gm + Gf12) - Vf * (Gm - Gf12) ) / ( (Gm + Gf12) + Vf * (Gm - Gf12) ));
G12_C_GPa(i) = G12 / 10^9;

end

% Graficar todos los resultados
figure;
plot(Vf_range, G12_A_GPa, '-o', 'LineWidth', 2, 'MarkerSize', 4,'DisplayName','Inciso A');
hold on;
plot(Vf_range, G12_B_GPa, '-*', 'LineWidth', 2, 'MarkerSize', 4,'DisplayName','Inciso B');
plot(Vf_range, G12_C_GPa, 'o--b', 'LineWidth', 2, 'MarkerSize', 4,'DisplayName','Inciso C');
xlabel('Vf (Fracción Volumétrica de la Fibra)');
ylabel('G12 (GPa)');
title('Gráfico de G12 vs Vf (Problema 4)');
grid on;
legend;
hold off;

%%
%Problema 5
Em=4.62e9;
Ef1=233e9;
Ef2=23.1e9;
Gf12=8.96e9;
vm=0.360;
vf12=0.200;
vf23=0.400;
Gf23=8.27e9;
Vf=0.6;
Vm=1-Vf;
am=41.4e-6;
af1=-0.540e-6;
af2=10.10*10^-6;
Gm=Em/(2*(1+vm));
E1 = (Ef1*Vf) + (Em*Vm);
E1_A_GPa = E1 / 10^9;
E2 = (Vf/Ef2) + (Vm/Em);
E2_A = (E2)^-1;
E2_A_GPa= E2_A / 10^9;
v12=(Vf*vf12)+(Vm*vm);
G12=(Vf/Gf12)+(Vm/Gm);
G12=(G12)^-1;
G12_A_GPa = G12 / 10^9;

Alpha1=((af1*Ef1*Vf)+(am*Em*Vm))/((Ef1*Vf)+(Em*Vm));
Alpha2=(Vf*af2)+(Vm*am);
Alpha2_1=(af2 - (Em/E1)*vf12*(am - af1)*Vm)*Vf + ...
         (am + (Ef1/E1)*vm*(am - af1)*Vf)*Vm;

% Imprimir los resultados en la ventana de comandos
fprintf('Resultados de la Actividad Fundamental 5 Problema 5\n');
fprintf('A) Alpha1: %.4e 1/K\n', Alpha1);
fprintf('B) Alpha2: %.4e 1/K (Regla de Mezclas)\n', Alpha2);
fprintf('C) Alpha2: %.4e 1/K\n', Alpha2_1);