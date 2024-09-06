close all
clear

%Young's Moduli for the Matrix 
YoungF=14.8;
YoungM=3.45;
vm = 0.36;

%Shear Moduli for Matrix
G12=4.40;
G13=4.40;
G23=3.20;

%Poisson Moduli for Matriz
v12=0.248;
v13=0.248;
v23=0.458;

Vf=1-vm;
 Vm=1-Vf;

Young1=(YoungM*Vm)+(YoungF*Vf)
Young2=(YoungM*YoungF)/((Vm*YoungF)+(Vf*YoungM))
Young3=Young2 %E2 == E3 fo this example

v21=(Young2*v12)/Young1;
v31=(Young3*v13)/Young1;
v32=(Young3*v23)/Young2;

%Matrix S
S=[1/Young1 -v21/Young2 -v31/Young3 0 0 0
    -v12/Young1 1/Young2 -v32/Young3 0 0 0
    -v13/Young1 -v23/Young2 1/Young3 0 0 0
    0 0 0 1/G23 0 0
    0 0 0 0 1/G13 0
    0 0 0 0 0 1/G12]

%Stresses Vector
sigma2=100/(60*60);

sigma=[0 sigma2 0 0 0 0]

sigma=sigma'

epsilon=S*sigma

d1=epsilon(1)*60
d2=epsilon(2)*60
d3=epsilon(3)*60