clc;

% Given material properties
V_f = 0.5;
V_m = 1 - V_f ;
E_f = 303.83e9; % GPa to Pa
E_m = 6.17e9; % GPa to Pa
nu12 = 0.248;
nu13 = 0.248;
nu23 = 0.458;
nu32 = nu23;
G23 = 3.20e9; % GPa to Pa
G12 = 4.40e9; % GPa to Pa
G13 = G12;

E1 = E_f * V_f + E_m * V_m

E2 = (V_f/E_f + V_m / E_m)^-1
E2 = G23 * 2 * (1 + nu23)
E2 = ((V_f/E_f + 0.5 * V_m) / (V_f + 0.5*V_m))^-1

E3 = nu32*E2/nu23

S = [
    1/E1, -nu12/E1, -nu13/E1, 0, 0, 0;
   -nu12/E2, 1/E2, -nu23/E2, 0, 0, 0;
   -nu13/E3, -nu23/E3, 1/E3, 0, 0, 0;
        0,      0,      0, 1/G23, 0, 0;
        0,      0,      0, 0, 1/G13, 0;
        0,      0,      0, 0, 0, 1/G12
]

% Stress in 2-direction
sigma2 = 100e3 / (60e-3 * 60e-3); % 100 kN, area 60x60 mm

strain = S * [0; sigma2; 0; 0; 0; 0];

