% Given material properties
V_f = 0.5;
E_f = 303.83e9; % GPa to Pa
E_m = 6.17e9; % GPa to Pa
nu12 = 0.248;
nu13 = 0.248;
nu23 = 0.458;
G23 = 3.20e9; % GPa to Pa
G12 = 4.40e9; % GPa to Pa
G13 = G12;

% Calculate effective Young's modulus
E1 = V_f * E_f + (1 - V_f) * E_m;
E2 = (E_f * E_m) / (V_f * E_m + (1 - V_f) * E_f);
E3 = E2; % Since E2 = E3 in the problem statement

% Stress in 2-direction
sigma2 = 100e3 / (60e-3 * 60e-3); % 100 kN, area 60x60 mm

% Compliance matrix
S = [
    1/E1, -nu12/E1, -nu13/E1, 0, 0, 0;
   -nu12/E2, 1/E2, -nu23/E2, 0, 0, 0;
   -nu13/E3, -nu23/E3, 1/E3, 0, 0, 0;
        0,      0,      0, 1/G23, 0, 0;
        0,      0,      0, 0, 1/G13, 0;
        0,      0,      0, 0, 0, 1/G12
];

% Strain vector
strain = S * [0; sigma2; 0; 0; 0; 0];

% Dimension changes
L0 = 60e-3; % Original length in meters
delta_L1 = L0 * strain(1);
delta_L2 = L0 * strain(2);
delta_L3 = L0 * strain(3);

% Display results
fprintf('Change in length in 1-direction: %.6f mm\n', delta_L1 * 1e3);
fprintf('Change in length in 2-direction: %.6f mm\n', delta_L2 * 1e3);
fprintf('Change in length in 3-direction: %.6f mm\n', delta_L3 * 1e3);
