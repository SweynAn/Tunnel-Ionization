function [ y ] = cNL( Z,ionizationPotential )
y=((2^(2*adjustedQuantum(Z,ionizationPotential)))/((adjustedQuantum(Z,ionizationPotential) ...
    *gamma(2*adjustedQuantum(Z,ionizationPotential))*gamma(1))));

end