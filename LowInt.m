function [ I ] = LowInt(ion,lambda,cycle, it)
% find the lowest possible laser intensity used in different gas
% cycle is how many cycles of ion you hope to use
% with given ion and wavelength
% mainly Argon, Helium and Neon
% Argon ==> Ar8+ 
% Neon ==> Ne8+
% Helium ==> He2+
% lambda ==> wavelength in nm
% it ==> iteration time you want to have, usually more than 1000 to
% ensure the value output is correct

% 2.42*10^-17 converts the value into atomic unit
% Using laser of wavelength lambda 
% omega = 2*pi*3*10^8/lambda
omega1=2*pi*3*10^8/lambda*10^9*2.42*10^(-17);
fs=10^(-15)/(2.42*10^(-17)); % atomic unit
Il0= 1*10^19 ; % Peak intensity Used in Argon initially   
tau = lambda*10^(-9)/(3*10^8)/(1.76*2.42*10^-17)*cycle;
Index = 0;
N=1;  % # of Ar
N1=0; % # of Ar+
N2=0; % # of Ar2+
N3=0; % # of Ar3+
N4=0; % # of Ar4+
N5=0;
N6=0;
N7=0;
N8=0;
N9=0;
N10=0;
N11=0;
N12=0;
N13=0;
N14=0;
N15=0;
N16=0;
N17=0;


adk=0;
derN=0;
derN1=0;
derN2=0;
derN3=0;
derN4=0;
derN5=0;
derN6=0;
derN7=0;
derN8=0;
derN9=0;
derN10=0;
derN11=0;
derN12=0;
derN13=0;
derN14=0;
derN15=0;
derN16=0;
derN17=0;

Cutoff=0;
%Switch the ion
switch ion(1:2)
    case 'Ar'
    while(Il0>0)
        % adkPeaks=omegaADK(E(tslow),1,.58,1,0); testing use
        % Used to Calculate ionization rate with certain population
        % N=1;N1=0;N2=0;N3=0;N4=0;adk=0;derN=0;derN1=0;derN2=0;derN3=0;derN4=0;
        
        % Assigned the value for N and omegaADK and cutoff Energy
        % Each step is pi/omega1/100*2.42 corresponding to half cycle in fs
        % 100fs is the total time range
        for i=1:ceil(100*fs/(pi/omega1/it*2.42))                        
        El0 = sqrt(Il0/(1*10^14)) * 0.053376 ;
        E= @(t) El0.*sech(t./tau).*abs(cos(omega1.*t));
        Inp = @(t) Il0.*(sech(t./tau)^2);
        %ceil rounds the elements of A to the nearest integer
    
        % Calculate out the value of Omega first
        
        adk(i)=omegaADK(E(2.42*pi/it/omega1*(i-1)-50*fs) ,1,.58,1,0);
        adk1(i)=omegaADK(E(2.42*pi/it/omega1*(i-1)-50*fs),2,1.01,1,0);
        adk2(i)=omegaADK(E(2.42*pi/it/omega1*(i-1)-50*fs),3,1.496,1,0);
        adk3(i)=omegaADK(E(2.42*pi/it/omega1*(i-1)-50*fs),4,2.1989,1,0);
        adk4(i)=omegaADK(E(2.42*pi/it/omega1*(i-1)-50*fs),5,2.758,1,0);
        adk5(i)=omegaADK(E(2.42*pi/it/omega1*(i-1)-50*fs),6,3.3459,1,0);
        adk6(i)=omegaADK(E(2.42*pi/it/omega1*(i-1)-50*fs),7,4.57069853,1,0);
        adk7(i)=omegaADK(E(2.42*pi/it/omega1*(i-1)-50*fs),8,5.27426471,1,0);
        adk8(i) =omegaADK(E(2.42*pi/it/omega1*(i-1)-50*fs),9, 15.53125,1,0);
        adk9(i) =omegaADK(E(2.42*pi/it/omega1*(i-1)-50*fs),10,17.5988971,1,0);
        adk10(i)=omegaADK(E(2.42*pi/it/omega1*(i-1)-50*fs),11,19.8147059,1,0);
        adk11(i)=omegaADK(E(2.42*pi/it/omega1*(i-1)-50*fs),12,22.7301471,1,0);
        adk12(i)=omegaADK(E(2.42*pi/it/omega1*(i-1)-50*fs),13,25.2242647,1,0);
        adk13(i)=omegaADK(E(2.42*pi/it/omega1*(i-1)-50*fs),14,27.7845588,1,0);
        adk14(i)=omegaADK(E(2.42*pi/it/omega1*(i-1)-50*fs),15,31.4253576,1,0);
        adk15(i)=omegaADK(E(2.42*pi/it/omega1*(i-1)-50*fs),16,33.7511029,1,0);
        adk16(i)=omegaADK(E(2.42*pi/it/omega1*(i-1)-50*fs),17,151.503151,1,0);
        adk17(i)=omegaADK(E(2.42*pi/it/omega1*(i-1)-50*fs),18,162.729029,1,0);  
        % Calculate the Change of population over time
        % [Might use other code/method to Calculate N]
        N(i+1) = N(i)-derN(i)*2.42*pi/omega1/it;
        % similarly in N1 ... N2 ... N3 ... etc
        % N(i+1)=N(i)*exp(2.42*(-adk(i))*pi/omega1/it);
        % Method to take the derivative
        derN(i+1)=adk(i)*N(i);
        derN1(i+1)=adk(i)*N(i) - adk1(i)*N1(i);
        derN2(i+1)=adk1(i)*N1(i)-adk2(i)*N2(i);
        derN3(i+1)=adk2(i)*N2(i)-adk3(i)*N3(i);
        derN4(i+1)=adk3(i)*N3(i)-adk4(i)*N4(i);
        derN5(i+1)=adk4(i)*N4(i)-adk5(i)*N5(i);
        derN6(i+1)=adk5(i)*N5(i)-adk6(i)*N6(i);
        derN7(i+1)=adk6(i)*N6(i)-adk7(i)*N7(i);
        % Update to higher ionization rate
        derN8(i+1)=adk7(i)*N7(i)-adk8(i)*N8(i);
        derN9(i+1)=adk8(i)*N7(i)-adk9(i)*N9(i);
        derN10(i+1)=adk9(i)*N9(i)-adk10(i)*N10(i);
        derN11(i+1)=adk10(i)*N10(i)-adk11(i)*N11(i);
        derN12(i+1)=adk11(i)*N11(i)-adk12(i)*N12(i);
        derN13(i+1)=adk12(i)*N12(i)-adk13(i)*N13(i);
        derN14(i+1)=adk13(i)*N13(i)-adk14(i)*N14(i);
        derN15(i+1)=adk14(i)*N14(i)-adk15(i)*N15(i);
        derN16(i+1)=adk15(i)*N15(i)-adk16(i)*N16(i);
        derN17(i+1)=adk16(i)*N16(i);
        % Method for using the 
        N1(i+1)=N1(i)+derN1(i)*2.42*pi/omega1/it;
        N2(i+1)=N2(i)+derN2(i)*2.42*pi/omega1/it;
        N3(i+1)=N3(i)+derN3(i)*2.42*pi/omega1/it;
        N4(i+1)=N4(i)+derN4(i)*2.42*pi/omega1/it;
        N5(i+1)=N5(i)+derN5(i)*2.42*pi/omega1/it;
        N6(i+1)=N6(i)+derN6(i)*2.42*pi/omega1/it;
        N7(i+1)=N7(i)+derN7(i)*2.42*pi/omega1/it;
        % Update to higher ion 
        N8(i+1)=N8(i)+derN8(i)*2.42*pi/omega1/it;
        N9(i+1)=N9(i)+derN9(i)*2.42*pi/omega1/it;
        N10(i+1)=N10(i)+derN10(i)*2.42*pi/omega1/it;
        N11(i+1)=N11(i)+derN11(i)*2.42*pi/omega1/it;
        N12(i+1)=N12(i)+derN12(i)*2.42*pi/omega1/it;
        N13(i+1)=N13(i)+derN13(i)*2.42*pi/omega1/it;
        N14(i+1)=N14(i)+derN14(i)*2.42*pi/omega1/it;
        N15(i+1)=N15(i)+derN15(i)*2.42*pi/omega1/it;
        N16(i+1)=N16(i)+derN16(i)*2.42*pi/omega1/it;
        N17(i+1)=N17(i)+derN17(i)*2.42*pi/omega1/it;
  
        end
        
        % another way is to calculate the index of the condition
        % get the least intensity for ion
        Index = ceil(length(N)/2) ;
        switch ion
            case 'Ar'
                if N1(Index)>0.1 && N(Index)<0.01
                   I = Il0;
                   break;
                else
                   Il0 = Il0*1.01;  
                end
            case 'Ar+'
                if N2(Index)>0.1 && N1(Index)<0.01
                   I = Il0;
                   break;
                else
                   Il0 = Il0*1.01;  
                end            
            case 'Ar2+'
                if N3(Index)>0.1 && N2(Index)<0.01
                   I = Il0;
                   break;
                else
                   Il0 = Il0*1.01;  
                end
            case 'Ar3+'
                if N4(Index)>0.1 && N3(Index)<0.01
                   I = Il0;
                   break;
                else
                   Il0 = Il0*1.01;  
                end
           % update to higher ions
            case 'Ar4+'
                if N5(Index)>0.1 && N4(Index)<0.01
                   I = Il0;
                   break;
                else
                   Il0 = Il0*1.01;  
                end
            case 'Ar5+'
                if N6(Index)>0.1 && N5(Index)<0.01
                   I = Il0;
                   break;
                else
                   Il0 = Il0*1.01;  
                end            
            case 'Ar6+'
                if N7(Index)>0.1 && N6(Index)<0.01
                   I = Il0;
                   break;
                else
                   Il0 = Il0*1.01;  
                end
            case 'Ar7+'
                if N8(Index)>0.1 && N7(Index)<0.01
                   I = Il0;
                   break;
                else
                   Il0 = Il0*1.01;  
                end      
        end
    end
    case 'He'     
            while(Il0>0)
            N=1;N1=0;N2=0;adk=0;derN=0;derN1=0;derN2=0;
            for i=1:ceil(100*fs/(pi/omega1/100*2.42)) 
                El0 = sqrt(Il0/(1*10^14)) * 0.053376;
                E= @(t) El0.*sech(t./tau).*abs(cos(omega1.*t));
                Inp = @(t) Il0.*(sech(t./tau)^2);
    
                % Calculate out the value of Omega first
                adk(i) = omegaADK(E(2.42*pi/100/...
                     omega1*(i-1)-50*fs),1,.904,1,0); % He --> He+
                adk1(i)= omegaADK(E(2.42*pi/100/...
                     omega1*(i-1)-50*fs),2,2,1,0); % He+ --> He2+
 
                % Calculate the Change of population over time
                % adk multiply the value of each step pi/omega1/100
                % N(i+1)=N(i)*exp(-adk(i)*2.42*pi/omega1/100); %???
                % Method to take the derivative to calculate 
                derN(i+1) =adk(i)*N(i);
                derN1(i+1)=adk(i)*N(i)-adk1(i)*N1(i);
                derN2(i+1)=adk1(i)*N1(i);
  
                % Method for using the
                N(i+1) = N(i)-derN(i)*2.42*pi/omega1/100;
                N1(i+1)=N1(i)+derN1(i)*2.42*pi/omega1/100;
                N2(i+1)=N2(i)+derN2(i)*2.42*pi/omega1/100;

                % Final rate should be the omegaADK
                % multiply the value of population
                RateN1(i+1)=adk1(i)*N1(i);
            end
            
            Index = ceil(length(N)/2);
            % To judge the break time
            switch ion
                case 'He'
                    if N2(Index)>0.1 && N1(Index)<0.01
                            I = Il0;
                            break;
                    else
                            Il0 = Il0*1.01;  
                    end                
                case 'He+'
                    if N2(Index)>0.1 && N1(Index)<0.01
                            I = Il0;
                            break;
                    else
                            Il0 = Il0*1.01;  
                    end
            end
            end
    case 'Ne'
    while(Il0>0)
        % adkPeaks=omegaADK(E(tslow),1,.58,1,0); testing use
        % Used to Calculate ionization rate with certain population
        N=1;N1=0;N2=0;N3=0;N4=0;adk=0;derN=0;derN1=0;derN2=0;derN3=0;derN4=0;
        
        % Assigned the value for N and omegaADK and cutoff Energy
        % Each step is pi/omega1/100*2.42 corresponding to half cycle in fs
        % 100fs is the total time range
        for i=1:ceil(100*fs/(pi/omega1/100*2.42))                        
        El0 = sqrt(Il0/(1*10^14)) * 0.053376 ;
        E= @(t) El0.*sech(t./tau).*abs(cos(omega1.*t));
        Inp = @(t) Il0.*(sech(t./tau)^2);
        %ceil rounds the elements of A to the nearest integer
    
        % Calculate out the value of Omega first
        adk(i) =omegaADK(E(2.42*pi/100/omega1*(i-1)-50*fs),1,.79281618,1,0);
        adk1(i)=omegaADK(E(2.42*pi/100/omega1*(i-1)-50*fs),2,1.506,1,0);
        adk2(i)=omegaADK(E(2.42*pi/100/omega1*(i-1)-50*fs),3,2.33272,1,0);
        adk3(i)=omegaADK(E(2.42*pi/100/omega1*(i-1)-50*fs),4,3.57,1,0);
        adk4(i)=omegaADK(E(2.42*pi/100/omega1*(i-1)-50*fs),5,4.64007,1,0);
  
        % Calculate the Change of population over time
        % [Might use other code/method to Calculate N]
        % [N(i+1) = N(i)-derN(i)*2.42*pi/omega1/100] 
        % similarly in N1 ... N2 ... N3 ... etc
        N(i+1)=N(i)*exp(2.42*(-adk(i))*pi/omega1/100);
  
        % Method to take the derivative
        derN(i+1)=adk(i)*N(i);
        derN1(i+1)=adk(i)*N(i)-adk1(i)*N1(i);
        derN2(i+1)=adk1(i)*N1(i)-adk2(i)*N2(i);
        derN3(i+1)=adk2(i)*N2(i)-adk3(i)*N3(i);
        derN4(i+1)=adk3(i)*N3(i);
  
        % Method for using the 
        N1(i+1)=N1(i)+derN1(i)*2.42*pi/omega1/100;
        N2(i+1)=N2(i)+derN2(i)*2.42*pi/omega1/100;
        N3(i+1)=N3(i)+derN3(i)*2.42*pi/omega1/100;
        N4(i+1)=N4(i)+derN4(i)*2.42*pi/omega1/100;
  
        % Final rate should be the omegaADK multiply the value of population
        RateN1(i+1)=adk1(i)*N1(i);
        RateN2(i+1)=adk2(i)*N2(i);
        RateN3(i+1)=adk3(i)*N3(i);
        
        end
        
        Index = ceil(length(N)/2) ;
        % another way is to calculate the index of the condition
        % get the least intensity for ion
        switch ion
            case 'Ne'
                if N1(Index)>0.1 && N(Index)<0.01
                   I = Il0;
                   break;
                else
                   Il0 = Il0*1.01;  
                end
            case 'Ne+'
                if N2(Index)>0.1 && N1(Index)<0.01
                   I = Il0;
                   break;
                else
                   Il0 = Il0*1.01;  
                end            
            case 'Ne2+'
                if N3(Index)>0.1 && N2(Index)<0.01
                   I = Il0;
                   break;
                else
                   Il0 = Il0*1.01;  
                end
            case 'Ne3+'
                if N4(Index)>0.1 && N3(Index)<0.01
                   I = Il0;
                   break;
                else
                   Il0 = Il0*1.01;  
                end
        end   
    end

end

