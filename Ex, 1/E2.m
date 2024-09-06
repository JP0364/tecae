% AF4_Book


function y = E2Modified(Vf,E2f,Em,Eta,NU12f,NU21f,NUm,E1f,p)
    %E2Modified Thisfunction returns Young’s modulus in the
    % transverse direction.Itsinputarenine values:
    % Vf- fiber volume fraction
    % E2f- transverse Young’s modulus of the fiber
    % Em- Young’s modulus of thematrix
    % Eta- stress-partitioning factor
    % NU12f- Poisson’s ratioNU12of thefiber
    % NU21f- Poisson’s ratioNU21of thefiber
    % NUm- Poisson’s ratioof thematrix
    % E1f- longitudinal Young’s modulus of the fiber
    % p- parameter usedto determinewhich equation to use:
    % p = 1- useequation (3.4)
    % p = 2- useequation (3.9)
    % p = 3- useequation (3.10)
    % p = 4- usethemodified formula using (3.23)
    % Usethevalue zeroforanyargument notneeded
    % in thecalculations.
    Vm=1-Vf;
    if p==1
    y = 1/(Vf/E2f+ Vm/Em);
    elseif p == 2;
    y = 1/((Vf/E2f+ Eta*Vm/Em)/(Vf + Eta*Vm));
    elseif p == 3;
    deno= E1f*Vf + Em*Vm;
    etaf= (E1f*Vf + ((1-NU12f*NU21f)*Em + NUm*NU21f...
    *E1f)*Vm)/deno;
    etam = (((1-NUm*NUm)*E1f- (1-NUm*NU12f)*Em)*Vf...
    + Em*Vm)/deno;
    y = 1/(etaf*Vf/E2f + etam*Vm/Em);
    elseif p == 4
    EmPrime = Em/(1- NUm*NUm);
 