% This function returns the ionization rate
% It takes input E for electric field amplitude
% Z : next ion number. For a neutral atom this is 1
% iP: ionization potential eV/27.2 in atomic unit
% l : angular quantum number
% m : magnetic quantum number

function [ y ] = omegaADK( E,Z,iP,l,m )
y=(3.*E./(pi.*(2.*iP).^(3/2))).^(1/2).*cNL(Z,iP) ...
    .*fLM(l,m).*iP.*((2.*(2.*iP) ...
    .^(3/2))./E).^((2.*(Z)/(2.*iP).^(1/2))-1-abs(m)) ...
    .*exp((-2.*(2.*iP).^(3/2))./(3.*E));
end