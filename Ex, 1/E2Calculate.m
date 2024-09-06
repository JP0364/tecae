function [y] = E2Calculate(Vf,Vm,Ef1,Ef2,Em,eta,NUf12,NUf21,NUm,method)
    % E2Calculate calculates the transverse modulus E2 using the selected
    % method.
    % Inputs:
    % Vf: Volume fraction of the fiber material
    % Vm: Volume fraction of the matrix material
    % Ef1: Young's modulus of the first material in Pascals (Pa)
    % Ef2: Young's modulus of the fiber material in Pascals (Pa)
    % Em: Young's modulus of the matrix material in Pascals (Pa)
    % eta: Poisson's ratio of the matrix material
    % NUf12: Poisson's ratio between the first and second materials
    % NUf21: Poisson's ratio between the second and first materials
    % method: Method to calculate the transverse modulus E2
        % 'simple': Simple rule-of-mixtures formula
        % 'modified': Modified rule-of-mixtures formula
        % 'alternative': Alternative rule-of-mixtures formula
    % Outputs:
    % y: Transverse modulus E2 in Pa

    % Volume fraction of the matrix material
    if Vm ~= 0
        Vf = 1 - Vm;
    elseif Vf ~= 0
        Vm = 1 - Vf;
    end

    % Calculate the transverse modulus E2 using the selected method
    switch method
        case 'simple'
            y = 1 ./ (Vf./Ef2 + Vm./Em); % Pa
        case 'modified'
            y = 1 ./ ((Vf./Ef2 + eta.*Vm./Em) ./ (Vf + eta.*Vm)); % Pa
        case 'alternative'
            nf = (Ef1.*Vf + ((1 - NUf12.*NUf21)*Em + NUm.*NUf21.*Ef1).*Vm) / (Ef1.*Vf + Em.*Vm);
            nu_f = (((1 - NUm.*NUm).*Ef1 - (1 - NUm.*NUf12).*Em).*Vf + Em.*Vm) / (Ef1.*Vf + Em.*Vm);

            y = 1 ./ ((nf .* Vf) ./ Ef2 + (nu_f .* Vm) ./ Em); % Pa
    end
end
