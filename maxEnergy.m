function [ y ] = maxEnergy(ion,fieldIntensity,drivingWavelength )
% Returns the maximum energy (in eV) photon that can be created with HHG given an
% ionization potential, field intensity, and drivng wavelength.
% Ionization potential in eV.
% Field intensity in W/cm^2.
% Driving wavelength in nm. 
% Noble gas ionization potential:
% https://en.wikipedia.org/wiki/Ionization_energies_of_the_elements_(data_page)
% and follow the cutoff rule for high harmonic generation: 
% https://en.wikipedia.org/wiki/High_harmonic_generation

switch ion
    
    % ionization potential for Argon        || /27.2
    case  'Ar' 
        ionizationPotential = 15.8 ;     %
    case  'Ar+'
        ionizationPotential = 27.6 ;
    case  'Ar2+'
        ionizationPotential = 40.7 ;
    case  'Ar3+'
        ionizationPotential = 59.81 ; 
    case  'Ar4+'
        ionizationPotential = 75.02 ;    % 2.758088
    case  'Ar5+'
        ionizationPotential = 91.009 ;    % 3.3459
    case  'Ar6+'
        ionizationPotential = 124.323 ;   % 4.57069853
    case  'Ar7+'
        ionizationPotential = 143.460 ;   % 5.27426471
    case  'Ar8+'
        ionizationPotential = 422.45;     % 15.53125
    case  'Ar9+'
        ionizationPotential = 478.69 ;    % 17.5988971
    case  'Ar10+'
        ionizationPotential = 538.96 ;    % 19.8147059
    case  'Ar11+'
        ionizationPotential = 618.26 ;    % 22.7301471
    case  'Ar12+'
        ionizationPotential = 686.10 ;    % 25.2242647
    case  'Ar13+'
        ionizationPotential = 755.74 ;    % 27.7845588
    case  'Ar14+'
        ionizationPotential = 854.77 ;    % 31.4253576
    case  'Ar15+'
        ionizationPotential = 918.03 ;    % 33.7511029
    case  'Ar16+'
        ionizationPotential = 4120.8857 ; % 151.503151
    case  'Ar17+'
        ionizationPotential = 4426.2296 ; % 162.729029
    % ionization potential for heilium
    case  'He' 
        ionizationPotential = 24.6 ;
    case  'He+'
        ionizationPotential = 54.4 ;
    % ionization potential for Neon
    case  'Ne'
        ionizationPotential = 21.5646 ;
    case  'Ne+' 
        ionizationPotential = 40.96328 ;
    case  'Ne2+'
        ionizationPotential = 63.45	;
    case  'Ne3+'
        ionizationPotential = 97.12 ;
    case  'Ne4+'
        ionizationPotential = 126.21 ;
    case  'Ne5+'
        ionizationPotential = 157.93 ;
    case  'Ne6+'
        ionizationPotential = 207.2759;
    case  'Ne7+'
        ionizationPotential = 239.0989;
end

% Final Calculation for Energy
y=ionizationPotential+3.17*9.337*10^(-20)*(fieldIntensity)*(drivingWavelength)^2;

end

