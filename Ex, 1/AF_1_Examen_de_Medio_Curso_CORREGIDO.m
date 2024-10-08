% %%                               AF 6                                    %%
clc

fprintf('\n');
disp('Materiales Compuestos 2')
disp('EMC');
disp('Jesus Antonio Ramirez Alpizar - 2077851');
disp('Fecha: 02 / Octubre / 2024');
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

%% Problem 5 %%
clear
Ef1 = 230E9; % Pa
Ef2 = 25E9; % Pa
nuf12 = 0.22;
Em = 3.3E9; % Pa
num = 0.35;
Vf = 0.53:0.01:0.57;
Vm = 1-Vf;

% The simple rule-of-mixtures formula:
% inv(E2) = Vf / Ef2 + Vm / Em
E2 = 1./(Vf./Ef2 + Vm./Em).*1e-09; % GPa

E1 = (Ef1.*Vf + Em.*Vm).*1e-09; % GPa


figure('Name', 'Problem 5');
subplot(2,1,1)
plot(Vf,E1,'b-', 'LineWidth', 1.5,'Marker', 'pentagram')
hold on
legend('E1')

subplot(2,1,2)
plot(Vf,E2,'b-', 'LineWidth', 1.5,'Marker', 'pentagram')
legend('E2')

disp('Problema 5')
fprintf('Para Vf = %.2f\n',Vf(2))
fprintf('E1 = %.4f GPa\n',E1(2))
fprintf('E2 = %.4f GPa\n\n\n\n',E2(2))


%% Problema 6
clear
Gf12 = 28.3E9; %GPa
Gm = 1.27E9; %GPa
Vf = 0.53:0.01:0.57;
Vm = 1-Vf;

% (a) The simple rule-of-mixtures formula:
% inv(G12) = Vf / Gf12 + Vm / Gm
G12a = 1 ./ (Vf./Gf12 + Vm./Gm).*1E-09; % GPa

% (b) The modified rule-of-mixtures formula with etaPrime = 0.6.
% inv(G12) = (Vf / Gf12 + etaPrime*Vm / Gm) / (Vf + etaPrime*Vm)
etaPrime = 0.6;
G12b = (1./((Vf./Gf12 + etaPrime.*Vm./Gm) ./ (Vf + etaPrime.*Vm))).*1E-9; % GPa

% (c) The elasticity formula.
G12c = Gm .* ((Gm + Gf12) - Vf.*(Gm - Gf12)) ./ ((Gm + Gf12) + Vf .* (Gm - Gf12)) .* 1E-09; % GPa

figure('Name', 'Problem 6');
plot(Vf,G12a,'b-', 'LineWidth', 1.5,'Marker', 'o')
hold on
plot(Vf,G12b,'k-', 'LineWidth', 1.5,'Marker', '+')
plot(Vf,G12c,'y-', 'LineWidth', 2,'Marker', 'x')
legend('ROM','Modified ROM','Elasticity Formula')
title 'Problem 6: Vf vs. G12 (GPa)'

selInd = length(Vf);
disp('Problema 6')
fprintf('Las condiciones se cumplen en las 3 aproximaciones cuando Vf = %.2f\n',Vf(selInd))
fprintf('ROM -> %.4f\n',G12a(selInd))
fprintf('Modified ROM -> %.4f\n',G12b(selInd))
fprintf('Elasticity formula -> %.4f\n\n\n\n\n',G12c(selInd))

%% Problema 7
clear
E1 = 181E09;   % Pa
E2 = 10.30E09; % Pa 
% E3 = E2;       % Transverse modulus E3 (Pa)
nu12 = 0.28;
G12 = 7.17E09; % Pa

sigmaX = 230E06;   % Pa
sigmaY = -50E06;   % Pa
tauXY =  -45E06;   % Pa
sigma = [sigmaX; sigmaY; tauXY]; % Stress tensor (Pa)

% nu21 = E2 / E1 * nu12; % Poisson's ratio nu21

thetaDeg = 0;                         % Orientation angles (Deg)
theta = deg2rad(thetaDeg);            % Orientation angles (Rad)

% Calculate the reduced compliance matrix
S = [1/E1     -nu12/E1  0; ...
    -nu12/E1  1/E2      0; ...
    0         0         1/G12]; % Compliance matrix (1/Pa)

m = cos(theta); % Sine of the orientation angle
n = sin(theta); % Cosine of the orientation angle
    
T = [m.^2,   n.^2,   2.*m.*n; 
         n.^2,   m.^2,  -2.*m.*n; 
        -m.*n,   m.*n,   m.^2-n.^2]; % Transformation matrix

S1 = T \ (S * T); % Transformed compliance matrix (1/Pa)
epsilon = S1*sigma; % Strain tensor (Pa)

% The strain in the 1-2 coordinate are:
disp('Problem 7')
disp('Results:');
disp('For the given stresses, the strain in the 1-2 cordinates are');
disp(['The strain in the x dimension of the element is ', num2str(epsilon(1, 1), '%.4f')]);
disp(['The strain in the y dimension of the element is ', num2str(epsilon(2, 1), '%.4f')]);
disp(['The GammaXY for the element is ', num2str(epsilon(3, 1), '%.4f')]);
fprintf('\n\n');

msg = sprintf('For the given stresses, the strain in the 1-2 cordinates are:\n');
msg = [msg, sprintf('The strain in the x dimension of the element is %.4f.\n',epsilon(1, 1))];
msg = [msg, sprintf('The strain in the y dimension of the element is %.4f.\n',epsilon(2, 1))];
msg = [msg, sprintf('The GammaXY for the element is %.4f MPa.\n\n\n\n', epsilon(3, 1))];
msgbox(msg, 'EMC - Problem 7');


%% Problema 8
clear
l = 0.5;    % Initial length (m)
sigmaX = -500E06:1E06:500E06; % Stress in the x-direction (Pa)
sigma = [sigmaX; zeros(length(sigmaX),1)'; zeros(length(sigmaX),1)']; % Stress tensor (Pa)
% (a) the fibers are aligned along the x-axis.
thetaA = 0; % Orientation angle 1 (deg)

thetaDeg = thetaA;  % Orientation angles (deg)
theta = deg2rad(thetaDeg);            % Orientation angles (rad)

% Material characteristics
E1 = 70E09;    % Longitudinal modulus E1 (Pa)
E2 = 4E09;     % Transverse modulus E2 (Pa)
% E3 = E2;       % Transverse modulus E3 (Pa)
nu12 = 0.25;  % Poisson's ratio nu12
% nu13 = nu12;   % Poisson's ratio nu13

% nu21 = E2 / E1 * nu12; % Poisson's ratio nu21
G12 = 1;

% Calculate the reduced compliance matrix
S = [1/E1     -nu12/E1  0; ...
    -nu12/E1  1/E2      0; ...
    0         0         1/G12]; % Compliance matrix (1/Pa)

m = cos(theta); % Sine of the orientation angle
n = sin(theta); % Cosine of the orientation angle
epsilon = zeros(3,length(sigmaX)); % Strain tensor (Pa) 

for i = 1:length(sigmaX)
    T = [m.^2,   n.^2,   2.*m.*n; 
         n.^2,   m.^2,  -2.*m.*n; 
        -m.*n,   m.*n,   m.^2-n.^2]; % Transformation matrix

    S1 = T \ (S * T); % Transformed compliance matrix (1/Pa)
    epsilon(:, i) = S1*sigma(:,i); % Strain tensor (Pa)
end

% Calculate the changes in the 50-mm dimensions of the element
deltaX = zeros(1,length(sigmaX));
deltaY = zeros(1,length(sigmaX));
deltaX(1, :) = (epsilon(1, :) * l) * 1E3; % Change in the x dimension (mm)
deltaY(2, :) = (epsilon(2, :) * l) * 1E3; % Change in the y dimension (mm)

if length(sigmaX) <= 6
    msg = cell(1, length(theta) * 5);
    msg{1} = sprintf('AF6 - Problem 1\n\n');
    for i = 1:length(sigmaX)
        disp(['For a stress of ', num2str(sigmaX(i).*1E-6, '%.1f'), ' MPa:']);
        disp(['The change in the x dimension of the element is ', num2str(deltaX(1, i), '%.4f'), ' mm.']);
        disp(['The change in the y dimension of the element is ', num2str(deltaY(2, i), '%.4f'), ' mm.']);
        disp(['The Gamma XY for the element is ', num2str(epsilon(3, i), '%.4f')]);
        fprintf('\n\n');
    
        msg{5 * (i - 1) + 2} = sprintf('For the orientation angle of %.1f MPa:', sigmaX(i)*1E-6);
        msg{5 * (i - 1) + 3} = sprintf('\nThe change in the x dimension of the element is %.4f mm.', deltaX(1, i));
        msg{5 * (i - 1) + 4} = sprintf('\nThe change in the y dimension of the element is %.4f mm.', deltaY(2, i));
        msg{5 * (i - 1) + 5} = sprintf('\nThe Gamma XY for the element is %.4f\n\n', epsilon(3, i));
    end
    msgbox(msg, 'Problem 8');
else
    
    % Plot changes in direction in function of angle
    figure('Name', 'Problem 8 - Change in Dimension vs. stress (sigma) [Assuming stress(y)=0]')
    % set(gcf, 'Position', get(0, 'Screensize'));
    plot(sigmaX*1E-6, deltaX(1, :), 'r', 'LineWidth', 2);
    hold on
    plot(sigmaX*1E-6, deltaY(2, :), 'b', 'LineWidth', 2);
    xlabel('Stress (sigma)');
    ylabel('Change in Dimension (mm)');
    title('Change in Dimension vs. Stress (sigma) [Assuming stress(y)=0]');
    legend('Change in X Dimension', 'Change in Y Dimension');
    grid on
end


%% Problema 8a
clear
L = 0.5;    % Initial length (m)

% (a) the fibers are aligned along the x-axis.
% thetaA = 0; % Orientation angle 1 (deg)

% thetaDeg = thetaA;  % Orientation angles (deg)
% theta = deg2rad(thetaDeg);            % Orientation angles (rad)

% Material characteristics
E1 = 70E09;    % Longitudinal modulus E1 (Pa)
E2 = 4E09;     % Transverse modulus E2 (Pa)
% E3 = E2;       % Transverse modulus E3 (Pa)
nu12 = 0.25;  % Poisson's ratio nu12
% nu13 = nu12;   % Poisson's ratio nu13

% nu21 = E2 / E1 * nu12; % Poisson's ratio nu21

G12 = 1;

% Calculate the reduced compliance matrix
S = [1/E1     -nu12/E1  0; ...
    -nu12/E1  1/E2      0; ...
    0         0         1/G12]; % Compliance matrix (1/Pa)

% m = cos(theta); % Sine of the orientation angle
% n = sin(theta); % Cosine of the orientation angle

Delta1 = 2.329e-3;
Delta2 = -.291e-3;

%Deformacion unitaria
epsilon1=Delta1/L;
epsilon2=(Delta2*2)/L;
Gamma12=0;
epsilon = [epsilon1 epsilon2 Gamma12];
% 
% T = [m.^2,   n.^2,   2.*m.*n; 
%      n.^2,   m.^2,  -2.*m.*n; 
%     -m.*n,   m.*n,   m.^2-n.^2]; % Transformation matrix
% 
% S1 = T \ (S * T); % Transformed compliance matrix (1/Pa)
sigma = epsilon / S;

disp('Problem 8')
disp('Results:');
    
disp('For the given deformation:');
disp(['The stress in the x dimension of the element is ', num2str(sigma(1, 1)*1e-6, '%.4f'), ' MPa.']);
disp(['The stress in the y dimension of the element is ', num2str(sigma(1, 2)*1e-6, '%.4f'), ' MPa.']);
disp(['The Tau for the element is ', num2str(sigma(1, 3), '%.4f') ' MPa.']);
fprintf('\n\n');

msg = sprintf('For the given deformation:\n');
msg = [msg, sprintf('The stress in the x dimension of the element is %.4f MPa.\n',sigma(1, 1)*1e-6)];
msg = [msg, sprintf('The stress in the y dimension of the element is %.4f MPa.\n',sigma(1, 2)*1e-6)];
msg = [msg, sprintf('The Tau for the element is %.4f MPa.\n\n\n\n', sigma(1, 3))];
msgbox(msg, 'EMC - Problem 8');


%% Problema 9
clear
l = 0.5;    % Initial length (m)
sigma = [230e6; 0; 0]; % MPa
% (a) the fibers are aligned along the x-axis.
thetaA = 20; % Orientation angle 1 (deg)

thetaDeg = thetaA;  % Orientation angles (deg)
theta = deg2rad(thetaDeg);            % Orientation angles (rad)

% Material characteristics
E1 = 63.6E09;    % Longitudinal modulus E1 (Pa)
E2 = 3.5E09;     % Transverse modulus E2 (Pa)
E3 = E2;       % Transverse modulus E3 (Pa)
nu12 = 0.254;  % Poisson's ratio nu12
nu13 = nu12;   % Poisson's ratio nu13

nu21 = E2 / E1 * nu12; % Poisson's ratio nu21

%Assumed shear modulus of matrix
Gm = 1.27E09;   % Shear modulus of the matrix material in Pascals (Pa)
% Assuming Vf = 0.55
Vf = 0.55;   % Volume fraction of the fiber material
Vm = 1 - Vf; % Volume fraction of the matrix material
G12 = 1 / (Vf/E2 + Vm/Gm); % Shear modulus G12 (Pa)

% Calculate the reduced compliance matrix
S = [1/E1     -nu12/E1  0; ...
    -nu12/E1  1/E2      0; ...
    0         0         1/G12]; % Compliance matrix (1/Pa)

m = cos(theta); % Sine of the orientation angle
n = sin(theta); % Cosine of the orientation angle


T = [m.^2,   n.^2,   2.*m.*n; 
     n.^2,   m.^2,  -2.*m.*n; 
    -m.*n,   m.*n,   m.^2-n.^2]; % Transformation matrix

S1 = T \ (S * T); % Transformed compliance matrix (1/Pa)
epsilon = S1*sigma; % Strain tensor (Pa)

deltaX = (epsilon(1, 1) * l + l) * 1e03; % Change in the x dimension (mm)
deltaY = (epsilon(2, 1) * l + l) * 1e03; % Change in the y dimension (mm)



disp('Problem 9')
disp('Results:');
disp('For Kevlar, assuming: Gm = 1.27Gpa, Ef = E2 Vf = 0.55 :');
disp(['The change in the x dimension of the element is ', num2str(deltaX, '%.4f'), ' mm.']);
disp(['The change in the y dimension of the element is ', num2str(deltaY, '%.4f'), ' mm.']);
disp(['The Gamma XY for the element is ', num2str(epsilon(3, 1), '%.4f')]);
fprintf('\n');