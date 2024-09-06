% %%                               AF 6                                    %%
clc

fprintf('\n');
disp('Materiales Compuestos 2')
disp('Actividad Fundamental 6');
disp('Jesus Antonio Ramirez Alpizar - 2077851');
disp('Fecha: 9 / Septiembre / 2024');
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
% Consider a plane element of size 50mm × 50mm made of glass-reinforced 
% polymer composite material whose elastic constants are below. 
% The element is subjected to a tensile stress sigma = 100MPa in the x-direction. 
% Use MATLAB to calculate the strains and the deformed dimensions of the 
% element in the following three cases:
l = 50E-03;    % Initial length (m)
sigmaX = 100E06; % Stress in the x-direction (Pa)
sigma = [sigmaX; 0; 0]; % Stress tensor (Pa)
% (a) the fibers are aligned along the x-axis.
thetaA = 0; % Orientation angle θ (°)
% (b) the fibers are inclined to the x-axis with an orientation angle theta = 45◦.
thetaB = 45; % Orientation angle 45 (°)
% (c) the fibers are inclined to the x-axis with an orientation angle theta = −45◦
thetaC = -45; % Orientation angle -45 (°)

thetaDeg = [thetaA, thetaB, thetaC];  % Orientation angles (°)
theta = deg2rad(thetaDeg);            % Orientation angles (°)

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

nu21 = E2 / E1 * nu12; % Poisson's ratio nu21

% Calculate the reduced compliance matrix
S = [1/E1     -nu12/E1  0; ...
    -nu12/E1  1/E2      0; ...
    0         0         1/G12]; % Compliance matrix (1/Pa)

m = cos(theta); % Sine of the orientation angle
n = sin(theta); % Cosine of the orientation angle
epsilon = zeros(3,length(theta)); % Strain tensor (Pa) 

for i = 1:length(theta)
    T = [m(i).^2,   n(i).^2,  2.*m(i).*n(i); 
          n(i).^2,   m(i).^2,  -2.*m(i).*n(i); 
          -m(i).*n(i),  m(i).*n(i),   m(i).^2-n(i).^2]; % Transformation matrix

    S1 = T \ (S * T); % Transformed compliance matrix (1/Pa)
    epsiloni = S1*sigma; % Strain tensor (Pa)
    epsilon(:, i) = epsiloni; % Strain tensor (Pa)
end

% Calculate the changes in the 50-mm dimensions of the element
deltaX = zeros(1,length(theta));
deltaY = zeros(1,length(theta));
deltaX(1, :) = (epsilon(1, :) * l + l) * 1e03; % Change in the x dimension (mm)
deltaY(2, :) = (epsilon(2, :) * l + l) * 1e03; % Change in the y dimension (mm)

disp('AF6 - Problem 1')
disp('Results:');
if length(theta) <= 6
    msg = cell(1, length(theta) * 5);
    msg{1} = sprintf('AF6 - Problem 1\n\n');
    for i = 1:length(theta)
        disp(['For the orientation angle of ', num2str(thetaDeg(i), '%.1f'), ' degrees:']);
        disp(['The change in the x dimension of the element is ', num2str(deltaX(1, i), '%.4f'), ' mm.']);
        disp(['The change in the y dimension of the element is ', num2str(deltaY(2, i), '%.4f'), ' mm.']);
        disp(['The Gamma XY for the element is ', num2str(epsilon(3, i), '%.4f')]);
        fprintf('\n\n');
    
        msg{5 * (i - 1) + 2} = sprintf('For the orientation angle of %.1f degrees:', thetaDeg(i));
        msg{5 * (i - 1) + 3} = sprintf('\nThe change in the x dimension of the element is %.4f mm.', deltaX(1, i));
        msg{5 * (i - 1) + 4} = sprintf('\nThe change in the y dimension of the element is %.4f mm.', deltaY(2, i));
        msg{5 * (i - 1) + 5} = sprintf('\nThe Gamma XY for the element is %.4f\n\n', epsilon(3, i));
    end
    msgbox(msg, 'AF6 - Problem 1');
else
    disp('Results plotted');
    % Removed for a shorter code
end

%% Problem 2 %%
% Consider a glass-reinforced  polymer composite lamina with the elastic 
% constants as given in Problem 2.7. Use MATLAB to plot the values of the 
% five elastic constants Ex, vxy, Ey, νyx, and Gxy as a function of the 
% orientation angle theta in the range −pi/2 ≤ θ ≤ pi/2.

% Orientation angle (rad)
theta = -pi/2:0.01:pi/2;         % Define theta by step size
% theta = linspace(-pi/2,pi/2,19); % Define theta by number of steps
thetaDeg = rad2deg(theta);

% Calculate the reduced compliance matrix
S = [1/E1     -nu12/E1  0; ...
    -nu12/E1  1/E2      0; ...
    0         0         1/G12]; % Compliance matrix (1/Pa)

m = cos(theta); % Sine of the orientation angle
n = sin(theta); % Cosine of the orientation angle

epsilon = zeros(3,length(theta)); % Strain tensor (Pa) 
S1 = zeros(3,3); % Transformed compliance matrix (1/Pa)

for i = 1:length(theta)
    T = [m(i).^2,   n(i).^2,  2.*m(i).*n(i); 
          n(i).^2,   m(i).^2,  -2.*m(i).*n(i); 
          -m(i).*n(i),  m(i).*n(i),   m(i).^2-n(i).^2]; % Transformation matrix

    S1 = T \ (S * T); % Transformed compliance matrix (1/Pa)
    epsiloni = S1*sigma; % Strain tensor (Pa)
    epsilon(:, i) = epsiloni; % Strain tensor (Pa)
end

% Calculate the elastic constants
Ex = E1 ./ (m.^4 + (E1 ./ G12 - 2.*nu12) .* n.^2.*m.^2 + E1./E2 .* n.^4);                       % Longitudinal modulus Ex (Pa)
vxy = nu12 .* (n.^4 + m.^4) - (1 + E1./E2 - E1./G12) .* n.^2.*m.^2;                             % Poisson's ratio nuXY
vxy = vxy ./ (m.^4 + (E1 ./ G12 - 2.*nu12) .* n.^2.*m.^2 + E1./E2 .* n.^2);                     % Poisson's ratio nuXY
Ey = E2 ./ (m.^4 + (E2 ./ G12 - 2.*nu21) .* n.^2.*m.^2 + E2./E1 .* n.^4);                       % Transverse modulus Ey (Pa)
vyx = nu12 .* (n.^4 + m.^4) - (1 + E2/E1 - E2/G12) .* n.^2.*m.^2;                               % Poisson's ratio nuYX
vyx = vyx ./ (m.^4 + (E2/G12 - 2*nu12) .* n.^2.*m.^2 + (E2./E1).* n.^2);                        % Poisson's ratio nuYX
Gxy = G12 ./ (n.^4 + m.^4 + 2.*((1 + 2.*nu12).*2.*G12 ./ E1 + 2.*G12./E2 - 1) .* n.^2.*m.^2);   % Shear modulus GXY (Pa)
elasticConstants = [Ex./1E09; vxy; Ey./1E09; vyx; Gxy./1E09];
elasticConstantsStr = {'Ex', 'vxy', 'Ey', 'vyx', 'Gxy'};
elasticConstantsUnits = {'GPa', '', 'GPa', '', 'GPa'};
elasticConstantsColors = {'r', 'b', 'g', 'm', 'k'};
elasticConstantsDescription = {'Longitudinal Modulus Ex', 'Poisson''s Ratio nuYX', 'Transverse Modulus Ey', 'Poisson''s Ratio nuYX', 'Shear Modulus GXY'};

% Calculate the changes in the 50-mm dimensions of the element
deltaX = zeros(1,length(theta));
deltaY = zeros(1,length(theta));
deltaX(1, :) = (epsilon(1, :) * l + l) * 1e03; % Change in the x dimension (mm)
deltaY(2, :) = (epsilon(2, :) * l + l) * 1e03; % Change in the y dimension (mm)

disp('AF6 - Problem 2')
disp('Results:');
if length(theta) <= 6
    disp('Display results')
    % Section removed for shorter code
else
    disp('Results plotted');

    % Plot changes in direction in function of angle
    figure('Name', 'Change in Dimension vs. Orientation Angle')
    set(gcf, 'Position', get(0, 'Screensize'));
    subplot(2, 1, 1)
    plot(thetaDeg, deltaX(1, :), 'r', 'LineWidth', 2);
    hold on
    plot(thetaDeg, deltaY(2, :), 'b', 'LineWidth', 2);
    xlabel('Orientation Angle (°)');
    ylabel('Change in Dimension (mm)');
    title('Change in Dimension vs. Orientation Angle');
    legend('Change in X Dimension', 'Change in Y Dimension');
    grid on

    % Plot Gamma XY in function of angle
    subplot(2, 1, 2)
    plot(thetaDeg, epsilon(3, :), 'g', 'LineWidth', 2);
    xlabel('Orientation Angle (°)');
    ylabel('Gamma XY');
    title('Gamma XY vs. Orientation Angle');
    grid on

    % Plot the elastic constants
    figure('Name', 'Elastic Constants')
    set(gcf, 'Position', get(0, 'Screensize'));
    for i = 1:5
        if i > 3
            subplot(2, 2, i-1)
        else
            subplot(2, 3, i)
        end
        plot(thetaDeg, elasticConstants(i, :), elasticConstantsColors{i}, 'LineWidth', 2);
        xlabel('Orientation Angle (°)');
        ylabel([elasticConstantsStr{i}, ' (', elasticConstantsUnits{i}, ')']);
        title([elasticConstantsDescription{i}, ' vs. Orientation Angle']);
        grid on
    end
end
