%%                               AF 4                                    %%
clc

fprintf('\n');
disp('Materiales Compuestos 2')
disp('Actividad Fundamental 4 -  Calculo de propiedades mecanicas');
disp('Jesus Antonio Ramirez Alpizar - 2077851');
disp('Fecha: 17/Agosto/2024');
fprintf('\n');

disp('Universidad Autonoma de Nuevo Leon');
disp('Facultad de Ingenieria Mecanica y Electrica');
disp('Instructor: Israel De Santiago Talavera')
disp('Ciclo: Ago-Dic 2024');
disp('Grupo: 001');
disp('Horario: Viernes V4 - V6');
fprintf('\n\n\n');
disp('---------------------------------------------------');
fprintf('\n\n\n');

close all

%% Problem 1 %%
% Consider a carbon/epoxy composite lamina with the following matrix and 
% fiber material properties
Ef2 = 14.8E09; % Young's modulus of the fiber material in Pascals (Pa)
Em = 3.45E09;  % Young's modulus of the matrix material in Pascals (Pa)
nu_m = 0.36;   % Poisson's ratio of the matrix material

% Use MATLAB to calculate the transverse modulus E2 using the following 
% three methods (use Vf = 0.65):
Vf = 0.65; % Volume fraction of the fiber material

%  (a) The simple rule-of-mixtures formula.
%  (b) The modified rule-of-mixtures formula with eta = 0.5.
%  (c) The alternative rule-of-mixtures formula. 
%      For this case, use Ef1 = 85.6GPa, νf12 = νf21 = 0.3.
 
Vm = 1 - Vf; % Volume fraction of the matrix material

% (a) Rules of mixtures formula: 
% inv(E2) = Vf / Ef2 + Vm / Em

disp('AF4 - Problem 1')
E2a = inv(Vf/Ef2 + Vm/Em); % Pa

% (b) Modified rule-of-mixtures formula with η = 0.5. 
% inv(E2) = (Vf/Ef2 + eta*Vm/Em) / (Vf + nVm)

eta = 0.5;
E2b = inv((Vf/Ef2 + eta*Vm/Em) / (Vf + eta*Vm)); % Pa

% (c) alternative rule-of-mixtures formula
% inv(E2) = (nf * Vf) / Ef2 + (nm * Vm) / Em

Ef1 = 85.6E09;  % Elastic modulus of the first material,
vf12 = 0.3;     % Poisson's ratio between the first and second materials,
vf21 = vf12;    % Poisson's ratio between the second and first materials,

nf = (Ef1*Vf + ((1- vf12*vf21)*Em + nu_m*vf21*Ef1)*Vm) / (Ef1*Vf + Em*Vm);
nu_f = (((1 - nu_m^2)*Ef1 - (1- nu_m*vf12)*Em)*Vf + Em*Vm) / (Ef1*Vf + Em*Vm);

E2c = inv((nf * Vf) / Ef2 + (nu_f * Vm) / Em); % Pa

E2 = [E2a; E2b; E2c] .* 1E-09; % GPa

disp('Results:');
disp(['(a) E2 using simple rule-of-mixtures formula: ', num2str(E2(1), '%.4f'), ' GPa']);
disp(['(b) E2 using modified rule-of-mixtures formula: ', num2str(E2(2), '%.4f'), ' GPa']);
disp(['(c) E2 using alternative rule-of-mixtures formula: ', num2str(E2(3), '%.4f'), ' GPa']);
fprintf('\n\n');

msg = sprintf('(a) E2 using simple rule-of-mixtures formula: %.4f GPa\n', E2(1));
msg = [msg, sprintf('(b) E2 using modified rule-of-mixtures formula: %.4f GPa\n', E2(2))];
msg = [msg, sprintf('(c) E2 using alternative rule-of-mixtures formula: %.4f GPa\n', E2(3))];
msgbox(msg, 'AF4 - Problem 1');


%% Problem 2 %%
% Consider the glass/epoxy composite lamina of Problem 1. Use MATLAB to 
% plot a graph of the transverse modulus E2 versus the fiber volume 
% fraction Vf for each one of the following cases. Use all values 
% of Vf ranging from 0 to 1 (in increments of 0.1).
Vf = 0:0.1:1;
Vm = 1 - Vf;

% (a) the simple rule-of-mixtures formula.
% (b) the modified rule-of-mixtures formula with eta = 0.4.
% (c) the modified rule-of-mixtures formula with eta = 0.5.
% (d) the modified rule-of-mixtures formula with eta = 0.6.
% (e) the alternative rule-of-mixtures formula with the values given in 
%     part (c) of Problem 1.
% Make sure that all five graphs appear on the same plot.

% (a) The simple rule-of-mixtures formula:
% inv(E2) = Vf / Ef2 + Vm / Em
E2a = (1 ./ (Vf ./ Ef2 + Vm ./ Em)) .* 1E-09; % GPa

% (b) the modified rule-of-mixtures formula with eta = 0.4.
eta = 0.4;
E2b = (1 ./ ((Vf ./ Ef2 + eta.*Vm./Em) ./ (Vf + eta.*Vm))) .* 1E-09; % GPa

% (c) the modified rule-of-mixtures formula with eta = 0.5.
eta = 0.5;
E2c = (1 ./ ((Vf ./ Ef2 + eta.*Vm./Em) ./ (Vf + eta.*Vm))) .* 1E-09; % GPa

% (d) the modified rule-of-mixtures formula with eta = 0.6.
eta = 0.6;
E2d = (1 ./ ((Vf ./ Ef2 + eta.*Vm./Em) ./ (Vf + eta.*Vm))) .* 1E-09; % GPa

% (e) the alternative rule-of-mixtures formula with the values given in 
%     part (c) of Problem 1.
% inv(E2) = (nf * Vf) / Ef2 + (nm * Vm) / Em
nf = (Ef1*Vf + ((1- vf12*vf21)*Em + nu_m*vf21*Ef1).*Vm) ./ (Ef1*Vf + Em.*Vm);
nu_f = (((1 - nu_m^2)*Ef1 - (1- nu_m*vf12)*Em).*Vf + Em.*Vm) ./ (Ef1*Vf + Em.*Vm);
E2e = (1 ./ ((nf .* Vf) ./ Ef2 + (nu_f .* Vm) ./ Em)) .* 1E-09; % GPa


disp('AF4 - Problem 2');
disp('Results plotted on ''Figure 1: Problem 2''');
fprintf('\n\n');
figure('Name', 'Problem 2');
plot(Vf, E2a, 'r-', 'LineWidth', 1.5);
hold on
plot(Vf, E2b, 'y--', 'LineWidth', 1.5);
plot(Vf, E2c, 'b:', 'LineWidth', 1.5);
plot(Vf, E2d, '-.', 'LineWidth', 1.5, 'Color', [0.5 0.5 1]);
plot(Vf, E2e, 'k-', 'LineWidth', 1.5, 'Marker', 'pentagram', 'MarkerSize', 4);

grid on
xlabel('Fiber Volume Fraction (Vf)');
ylabel('Transverse Modulus (E2) [GPa]');
title('Transverse Modulus (E2) vs. Fiber Volume Fraction (Vf)');
legend('a) Simple Rule of mixtures', 'b) Modified Rule of mixtures (eta = 0.4)', ...
    'c) Modified Rule of mixtures (eta = 0.5)', 'd) Modified Rule of mixtures (eta = 0.6)', ...
    'e) Alternative Rule of mixtures');

%% Problem 3 %%
% Consider a carbon/epoxy composite lamina with the following matrix and 
% fiber material properties
Gf12 = 28.3E09; % Shear modulus of the fiber material in Pascals (Pa)
Gm = 1.27E09;   % Shear modulus of the matrix material in Pascals (Pa)

% Use MATLAB to calculate the shear modulus G12 using the following 
% three methods (use Vf = 0.55):
Vf = 0.55;   % Volume fraction of the fiber material
Vm = 1 - Vf; % Volume fraction of the matrix material
% (a) the simple rule-of-mixtures formula.
% (b) the modified rule-of-mixtures formula with etaPrime = 0.6.
% (c) the elasticity formula.

% (a) The simple rule-of-mixtures formula:
% inv(G12) = Vf / Gf12 + Vm / Gm
G12a = 1 / (Vf/Gf12 + Vm/Gm)*1E-09; % GPa

% (b) The modified rule-of-mixtures formula with etaPrime = 0.6.
% inv(G12) = (Vf / Gf12 + etaPrime*Vm / Gm) / (Vf + etaPrime*Vm)
etaPrime = 0.6;
G12b = (inv((Vf/Gf12 + etaPrime*Vm/Gm) / (Vf + etaPrime*Vm)))*1E-09; % GPa

% (c) The elasticity formula.
G12c = Gm * ((Gm + Gf12) - Vf*(Gm - Gf12)) / ((Gm + Gf12) + Vf * (Gm - Gf12)) * 1E-09; % GPa

G12 = [G12a; G12b; G12c];

disp('AF4 - Problem 3');
disp('Results:');
disp(['(a) G12 using simple rule-of-mixtures formula: ', num2str(G12(1), '%.4f'), ' GPa']);
disp(['(b) G12 using modified rule-of-mixtures formula: ', num2str(G12(2), '%.4f'), ' GPa']);
disp(['(c) G12 using elasticity formula: ', num2str(G12(3), '%.4f'), ' GPa']);
fprintf('\n\n');

msg = sprintf('(a) G12 using simple rule-of-mixtures formula: %.4f GPa\n', G12(1));
msg = [msg, sprintf('(b) G12 using modified rule-of-mixtures formula: %.4f GPa\n', G12(2))];
msg = [msg, sprintf('(c) G12 using elasticity formula: %.4f GPa\n', G12(3))];
msgbox(msg, 'AF4 - Problem 3');

%% Problem 4 %%
% Consider the glass/epoxy composite lamina of Problem 3. 
% Use MATLAB to plot a graph of the shear modulus G12 versus the fiber 
% volume fraction Vf for each one of the following cases. 
% Use all values of Vf ranging from 0 to 1 (in increments of 0.1).
Vf = 0:0.1:1; % Volume fraction of the fiber material
Vm = 1 - Vf;  % Volume fraction of the matrix material

% (a) the simple rule-of-mixtures formula.
% (b) the modified rule-of-mixtures formula with etaPrime = 0.6.
% (c) the elasticity formula.
% Make sure that all three graphs appear on the same plot

% (a) The simple rule-of-mixtures formula:
% inv(G12) = Vf / Gf12 + Vm / Gm
G12a = (1 ./ (Vf ./ Gf12 + Vm ./ Gm)) .* 1E-09; % GPa

% (b) The modified rule-of-mixtures formula with etaPrime = 0.6.
% inv(G12) = (Vf / Gf12 + etaPrime*Vm / Gm) / (Vf + etaPrime*Vm)
etaPrime = 0.6;
G12b = (1 ./ ((Vf ./ Gf12 + etaPrime.*Vm ./ Gm) ./ (Vf + etaPrime.*Vm))) .* 1E-09; % GPa

% (c) The elasticity formula.
G12c = (Gm .* ((Gm + Gf12) - Vf.*(Gm - Gf12)) ./ ((Gm + Gf12) + Vf.*(Gm - Gf12))) .* 1E-09; % GPa

disp('AF4 - Problem 4');
disp('Results plotted on ''Figure 2: Problem 4''');
fprintf('\n\n');
figure('Name', 'Problem 4');

plot(Vf, G12a, 'r-', 'LineWidth', 1.5);
hold on
plot(Vf, G12b, 'b--', 'LineWidth', 1.5);
plot(Vf, G12c, 'k:', 'LineWidth', 1.5);

grid on
xlabel('Fiber Volume Fraction (Vf)');
ylabel('Shear Modulus (G12) [GPa]');
title('Shear Modulus (G12) vs. Fiber Volume Fraction (Vf)');
legend('a) Simple Rule of mixtures', 'b) Modified Rule of mixtures (etaPrime = 0.6)', ...
    'c) Elasticity Formula');

%% Problem 5 %%
% Consider the graphite-reinforced polymer composite lamina of Example 2.
Em = 4.62E09;   % Young's modulus of the matrix material in Pascals (Pa)
Ef1 = 233E09;   % Young's modulus of the first material in Pascals (Pa)
Ef2 = 23.1E09;  % Young's modulus of the second material in Pascals (Pa)
Gf12 = 8.96E09; % Shear modulus of the fiber material in Pascals (Pa)
eta_m = 0.36;   % Poisson's ratio of the matrix material
eta_f12 = 0.2;  % Poisson's ratio between the first and second materials
eta_f23 = 0.4;  % Poisson's ratio between the second and third materials
Gf23 = 8.27E09; % Shear modulus of the third material in Pascals (Pa)

Vf = 0.6; % Volume fraction of the fiber material
Vm = 1 - Vf; % Volume fraction of the matrix material

% Let the coefficients of thermal expansion for the matrix and fibers 
% be given as follows:
alpha_m = 41.4E-06; % Coefficient of thermal expansion of the matrix material (1/K)
alpha_f1 = -0.54E-06; % Coefficient of thermal expansion of the first material (1/K)
alpha_f2 = 10.10E-06; % Coefficient of thermal expansion of the second material (1/K)

% Use MATLAB to calculate alpha1 and alpha2 for the lamina. When calculating α2, 
% use the two formulas.

% Thermal expansion coefficient alpha1:
alpha1 = ((alpha_f1*Ef1*Vf) + (alpha_m*Em*Vm)) / (Ef1*Vf + Em*Vm); % 1/K

% Longitudinal modulus E1:
E1 = (Ef1*Vf) + (Em*Vm); % Pa

% Thermal expansion coefficient alpha2 using simple rule-of-mixtures relation:
alpha2a = alpha_f2 * Vf + alpha_m * Vm; % 1/K

% Thermal expansion coefficient alpha2:
alpha2b = Vf * (alpha_f2 - (Em / E1) * eta_f12 * (alpha_m - alpha_f1) * Vm) + ...
          Vm * (alpha_m + (Ef1 / E1) * eta_m * (alpha_m - alpha_f1) * Vf); % 1/K

disp('AF4 - Problem 5');
disp('Results:');
disp(['Alpha1: ', num2str(alpha1, '%.4e'), ' 1/K']);
disp(['Alpha2: ', num2str(alpha2a, '%.4e'), ' 1/K (Simple Rule of Mixtures)']);
disp(['Alpha2: ', num2str(alpha2b, '%.4e'), ' 1/K']);
fprintf('\n\n');

msg = sprintf('Alpha1: %.4e 1/K\n', alpha1);
msg = [msg, sprintf('Alpha2: %.4e 1/K  (Simple Rule of Mixtures)\n', alpha2a)];
msg = [msg, sprintf('Alpha2: %.4e 1/K\n', alpha2b)];
msgbox(msg, 'AF4 - Problem 5');

clear
