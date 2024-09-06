% %%                               AF 5                                    %%
% clc
% 
% fprintf('\n');
% disp('Materiales Compuestos 2')
% disp('Actividad Fundamental 5');
% disp('Jesus Antonio Ramirez Alpizar - 2077851');
% disp('Fecha: 24 / Agosto / 2024');
% fprintf('\n');    
% 
% disp('Universidad Autonoma de Nuevo Leon');
% disp('Facultad de Ingenieria Mecanica y Electrica');
% disp('Instructor: Israel De Santiago Talavera')
% disp('Ciclo: Ago-Dic 2024');
% disp('Grupo: 001');
% disp('Horario: Viernes V4 - V6');
% fprintf('\n\n\n');
% disp('---------------------------------------------------');
% fprintf('\n\n\n');
% 
% close all
% 
% %% Problem 1 %%
% % Consider a 40-mm cube made of glass-reinforced polymer composite 
% % material that is subjected to a compressive force of 150 kN perpendicular 
% % to the fiber direction, directed along the 3-direction. The cube is free 
% % to expand or contract. Use MATLAB to determine the changes in the 40-mm 
% % dimensions of the cube. The material constants for glass-reinforced 
% % polymer composite material are given as follows: 
% % A=-0.0191, -0.1056, 0.2467
% F = 150E03;    % Force (N)
% l = 40E-03;    % Initial length (m)
% 
% E1 = 50.0E09;  % Longitudinal modulus E1 (Pa)
% E2 = 15.2E09;  % Transverse modulus E2 (Pa)
% E3 = E2;       % Transverse modulus E3 (Pa)
% nu23 = 0.428;  % Poisson's ratio nu23
% nu12 = 0.254;  % Poisson's ratio nu12
% nu13 = nu12;   % Poisson's ratio nu13
% G23 = 3.28E09; % Shear modulus G23 (Pa)
% G12 = 4.70E09; % Shear modulus G12 (Pa)
% G13 = G12;     % Shear modulus G13 (Pa)
% 
% % Calculate the stress tensor
% sigma3 = F / (l^2); % Stress in the 3-direction (Pa)
% sigma = [0; 0; sigma3; 0; 0; 0]; % Stress tensor (Pa)
% 
% % Calculate the strain tensor
% S = [1/E1     -nu12/E1  -nu13/E1  0      0      0; ...
%     -nu12/E1  1/E2      -nu23/E2  0      0      0; ...
%     -nu13/E1  -nu23/E2  1/E3      0      0      0; ...
%     0         0         0         1/G23  0      0; ...
%     0         0         0         0      1/G13  0; ...
%     0         0         0         0      0      1/G12]; % Compliance matrix (1/Pa)
% 
% % Calculate the strain tensor
% epsilon = S * sigma; % Strain tensor (Pa)
% 
% % Calculate the changes in the 40-mm dimensions of the cube
% delta_l = epsilon(1) * l; % Change in the 40-mm dimension (m)
% delta_2 = epsilon(2) * l; % Change in the 40-mm dimension (m)
% delta_1 = epsilon(3) * l; % Change in the 40-mm dimension (m)
% 
% disp('AF5 - Problem 1')
% disp('Results:');
% disp(['The change in the first dimension of the cube is  ', num2str(delta_l*1000, '%.4f'), ' mm.']);
% disp(['The change in the second dimension of the cube is ', num2str(delta_2*1000, '%.4f'), ' mm.']);
% disp(['The change in the third dimension of the cube is  ', num2str(delta_1*1000, '%.4f'), ' mm.']);
% fprintf('\n\n');
% 
% msg = sprintf('The change in the first dimension of the cube is  %.4f mm.\n', delta_l*1000);
% msg = [msg, sprintf('The change in the second dimension of the cube is %.4f mm.\n', delta_2*1000)];
% msg = [msg, sprintf('The change in the third dimension of the cube is  %.4f mm.\n', delta_1*1000)];
% msgbox(msg, 'AF5 - Problem 1');
% 
% %% Problem 2 %%
% % Repeat Problem 2.7 if the cube is made of aluminum instead of 
% % glass-reinforced polymer composite material. The material constants  
% % for aluminum are E = 72.4GPa and ν = 0.300. Use MATLAB.
% E = 72.4E09; % Young's modulus E (Pa)
% nu = 0.300;  % Poisson's ratio nu
% 
% % Calculate the stress tensor
% sigma = [0; 0; sigma3; 0; 0; 0]; % Stress tensor (Pa)
% 
% % Calculate the strain tensor
% S = [1/E     -nu/E  -nu/E  0      0      0; ...
%     -nu/E  1/E     -nu/E  0      0      0; ...
%     -nu/E  -nu/E  1/E     0      0      0; ...
%     0      0      0      1/G23  0      0; ...
%     0      0      0      0      1/G13  0; ...
%     0      0      0      0      0      1/G12]; % Compliance matrix (1/Pa)
% 
% % Calculate the strain tensor
% epsilon = S * sigma; % Strain tensor (Pa)
% 
% % Calculate the changes in the 40-mm dimensions of the cube
% delta_l = epsilon(1) * l; % Change in the 40-mm dimension (m)
% delta_2 = epsilon(2) * l; % Change in the 40-mm dimension (m)
% delta_1 = epsilon(3) * l; % Change in the 40-mm dimension (m)
% 
% disp('AF5 - Problem 2')
% disp('Results:');
% disp(['The change in the first dimension of the cube is  ', num2str(delta_l*1000, '%.4f'), ' mm.']);
% disp(['The change in the second dimension of the cube is ', num2str(delta_2*1000, '%.4f'), ' mm.']);
% disp(['The change in the third dimension of the cube is  ', num2str(delta_1*1000, '%.4f'), ' mm.']);
% fprintf('\n\n');
% 
% msg = sprintf('The change in the first dimension of the cube is  %.4f mm.\n', delta_l*1000);
% msg = [msg, sprintf('The change in the second dimension of the cube is %.4f mm.\n', delta_2*1000)];
% msg = [msg, sprintf('The change in the third dimension of the cube is  %.4f mm.\n', delta_1*1000)];
% msgbox(msg, 'AF5 - Problem 2');
% 
% %% Example 2 %%
% % Consider a 60-mm cube made of graphite-reinforced polymer composite 
% % material that is subjected to a tensile force of 100 kN perpendicular 
% % to the fiber direction, directed along the 2-direction. The cube is 
% % free to expand or contract. Use MATLAB to determine the changes in 
% % the 60-mm dimensions of the cube. The material constants for 
% % graphite-reinforced polymer composite material are given as follows:
% l = 60E-03;    % Initial length (m)
% F = 100E03;    % Force (N)
% 
% % Material Constants
% Vf = 0.50;       % Volume fraction of the fiber material
% Vm = 1 - Vf;     % Volume fraction of the matrix material
% Ef = 303.83E09;  % Young's modulus of the first material (Pa)
% Em = 6.17E09;    % Young's modulus of the matrix material (Pa)
% nu23 = 0.458;    % Poisson's ratio nu23
% nu12 = 0.248;    % Poisson's ratio nu12
% nu13 = nu12;     % Poisson's ratio nu13
% G23 = 3.20E09;   % Shear modulus G23 (Pa)
% G12 = 4.40E09;   % Shear modulus G12 (Pa)
% G13 = G12;       % Shear modulus G13 (Pa)
% 
% % Calculate the stress tensor
% E1 = Em * Vm + Ef * Vf; % Longitudinal modulus E1 (Pa)
% E2 = ((Em *Ef) / ((Vm * Ef) + Vf * Em)); % Transverse modulus E2 (Pa)
% E3 = E2; % Transverse modulus E3 (Pa)
% nu21 = (E2 * nu12) / E1; % Poisson's ratio nu21
% nu31 = (E2 * nu13) / E1; % Poisson's ratio nu31
% nu32 = (E2 * nu23) / E2; % Poisson's ratio nu32
% 
% % Problem 3 %
% % When a fiber-reinforced composite material is heated or cooled, the 
% % material expands or contracts just like an isotropic material. This is 
% % deformation that takes place independently of any applied load. 
% % Let deltaT be the change in temperature and let alpha1, alpha2, and 
% % alpha3 be the coefficients of thermal expansion for the composite 
% % material in the 1, 2, and 3-directions, respectively. In this case, 
% % the stress-strain relation of and becomes as follows:
% 
% % Consider now the cube of graphite-reinforced polymer composite material 
% % of Example 2. Suppose the cube is heated 30ºC above some reference state. 
% % Given alpha1 = −0.01800 × 10−6/◦C and alpha2 = alpha3 = 24.3×10−6/◦C, 
% % use MATLAB to determine the changes in length of the cube in each one 
% % of the three directions
% deltaT = 30;           % Change in temperature (°C)
% alpha1 = -0.01800E-06; % Coefficient of thermal expansion alpha1 (1/°C)
% alpha2 = 24.3E-06;     % Coefficient of thermal expansion alpha2 (1/°C)
% alpha3 = alpha2;       % Coefficient of thermal expansion alpha3 (1/°C)
% 
% % Calculate the stress tensor
% S = [1/E1     -nu21/E2  -nu31/E3  0      0      0; ...
%     -nu21/E2  1/E2      -nu32/E3  0      0      0; ...
%     -nu31/E3  -nu32/E3  1/E3      0      0      0; ...
%     0         0         0         1/G23  0      0; ...
%     0         0         0         0      1/G13  0; ...
%     0         0         0         0      0      1/G12]; % Compliance matrix (1/Pa)
% 
% % Calculate the strain tensor
% sigma2 = F / (l^2); % Stress in the 2-direction (Pa)
% sigma = [0; sigma2; 0; 0; 0; 0]; % Stress tensor (Pa)
% 
% % Calculate the strain tensor
% epsilonE = S * sigma; % Strain tensor (Pa)
% epsilonT = [alpha1 * deltaT; alpha2 * deltaT; alpha3 * deltaT; 0; 0; 0]; % Strain tensor (Pa)
% 
% % Calculate the changes in the 60-mm dimensions of the cube
% delta_l = epsilonE(1) * l + epsilonT(1) * l; % Change in the 60-mm dimension (m)
% delta_2 = epsilonE(2) * l + epsilonT(2) * l; % Change in the 60-mm dimension (m)
% delta_1 = epsilonE(3) * l + epsilonT(3) * l; % Change in the 60-mm dimension (m)
% 
% disp('AF5 - Problem 3')
% disp('Results:');
% disp(['The change in the first dimension of the cube is  ', num2str(delta_l*1000, '%.4f'), ' mm.']);
% disp(['The change in the second dimension of the cube is ', num2str(delta_2*1000, '%.4f'), ' mm.']);
% disp(['The change in the third dimension of the cube is  ', num2str(delta_1*1000, '%.4f'), ' mm.']);
% fprintf('\n\n');
% 
% msg = sprintf('The change in the first dimension of the cube is  %.4f mm.\n', delta_l*1000);
% msg = [msg, sprintf('The change in the second dimension of the cube is %.4f mm.\n', delta_2*1000)];
% msg = [msg, sprintf('The change in the third dimension of the cube is  %.4f mm.\n', delta_1*1000)];
% msgbox(msg, 'AF5 - Problem 3');

clear all
% %% Problem 4 %%
% % Consider a plane element of size 50mm × 50mm made of glass-reinforced 
% % polymer composite material whose elastic constants are below. 
% % The element is subjected to a tensile stress sigma = 100MPa in the x-direction. 
% % Use MATLAB to calculate the strains and the deformed dimensions of the 
% % element in the following three cases:
l = 50E-03;    % Initial length (m)
sigmaX = 100E06; % Stress in the x-direction (Pa)
sigma = [sigmaX; 0; 0]; % Stress tensor (Pa)
% (a) the fibers are aligned along the x-axis.
thetaA = 0; % Orientation angle θ (°)
% (b) the fibers are inclined to the x-axis with an orientation angle theta = 45◦.
thetaB = 45; % Orientation angle 45 (°)
% (c) the fibers are inclined to the x-axis with an orientation angle theta = −45◦
thetaC = -45; % Orientation angle -45 (°)
theta = [thetaA, thetaB, thetaC]; % Orientation angles (°)
% theta = -180:.01:180;

% Material characteristics
E1 = 50.0E09;  % Longitudinal modulus E1 (Pa)
E2 = 15.2E09;  % Transverse modulus E2 (Pa)
E3 = E2;       % Transverse modulus E3 (Pa)
nu23 = 0.428;  % Poisson's ratio nu23
nu12 = 0.254;  % Poisson's ratio nu12
nu13 = nu12;   % Poisson's ratio nu13
G23 = 3.28E09; % Shear modulus G23 (Pa)
G12 = 4.70E09; % Shear modulus G12 (Pa)
G13 = G12;     % Shear modulus G13 (Pa)

% Calculate the reduced compliance matrix
S = [1/E1     -nu12/E1  0; ...
    -nu12/E1  1/E2      0; ...
    0         0         1/G12]; % Compliance matrix (1/Pa)

m = cos(theta.*pi./180); % Cosine of the orientation angle
n = sin(theta.*pi./180);% Sine of the orientation angle

epsilon = zeros(3,length(theta)); % Strain tensor (Pa) 

for i = 1:length(theta)
    T = [m(i).^2,   n(i).^2,  2.*m(i).*n(i); 
          n(i).^2,   m(i).^2,  -2.*m(i).*n(i); 
          -m(i).*n(i),  m(i).*n(i),   m(i).^2-n(i).^2]; % Transformation matrix

    S1 = T \ (S * T); % Transformed compliance matrix (1/Pa)
    epsiloni = S1*sigma; % Strain tensor (Pa)
    epsilon(:, i) = epsiloni; % Strain tensor (Pa)
    theta(i)

end

% Calculate the changes in the 50-mm dimensions of the element
deltaX = zeros(1,length(theta));
deltaY = zeros(1,length(theta));
deltaX(1, :) = (epsilon(1, :) * l + l) * 1e03; % Change in the x dimension (mm)
deltaY(2, :) = (epsilon(2, :) * l + l) * 1e03; % Change in the y dimension (mm)

disp('AF5 - Problem 4')
disp('Results:');
for i = 1:length(theta)
    disp(['For the orientation angle of ', num2str(theta(i)), ' degrees:']);
    disp(['The change in the x dimension of the element is ', num2str(deltaX(1, i), '%.4f'), ' mm.']);
    disp(['The change in the y dimension of the element is ', num2str(deltaY(2, i), '%.4f'), ' mm.']);
    disp(['The Gamma XY for the element is ', num2str(epsilon(3, i), '%.4f')]);
    fprintf('\n');
end

% % Plot changes in direction in function of angle
% figure
% subplot(2, 1, 1)
% plot(theta, deltaX(1, :), 'r', 'LineWidth', 2);
% hold on
% plot(theta, deltaY(2, :), 'b', 'LineWidth', 2);
% xlabel('Orientation Angle (°)');
% ylabel('Change in Dimension (mm)');
% title('Change in Dimension vs. Orientation Angle');
% legend('Change in X Dimension', 'Change in Y Dimension');
% grid on

% % Plot Gamma XY in function of angle
% subplot(2, 1, 2)
% plot(theta, epsilon(3, :), 'g', 'LineWidth', 2);
% xlabel('Orientation Angle (°)');
% ylabel('Gamma XY');
% title('Gamma XY vs. Orientation Angle');
