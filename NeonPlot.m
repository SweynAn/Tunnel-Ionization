
% ==================================================
% Reproduce the data in PRL 103, 143901 (2009)
% Calculation of ionization Rate, ionizatio fraction
% with Cutoff Energy
% ==================================================

clc; clear; % Clean all the data in previous section
fs=10^(-15)/(2.42*10^(-17)); % atomic unit
Il0=1.8952e+16 ; % Peak intensity Used in Argon
El0=IE(Il0); % Amplitude for the electric field 
          % Calculated for the laser intensity
          % Intensity of 1* 10^14 W/cm^2


% 2.42*10^-17 converts the value into atomic unit
% Using laser of wavelength 800nm 
lambda = 343;
omega1=2*pi*3*10^8/lambda*10^9*2.42*10^(-17);


% 19fs here used in full pulse 
% tslow calculate the half cycle value
% 10.8 Sech^2 shaped pulses gives that 
% The full width at half-maximum pulse 
% duration is approx 1.76 times the parameter
% Inp is the laser intensity profile 
tau = lambda*10^(-9)/(3*10^8)/(1.76*2.42*10^-17)*7;
E= @(t) El0.*sech(t./tau).*abs(cos(omega1.*t));
Inp = @(t) Il0.*(sech(t./tau)^2); 
tslow=-50*fs:pi/omega1/100*2.42:50*fs;

% adkPeaks=omegaADK(E(tslow),1,.58,1,0); testing use
% Used to Calculate ionization rate with certain population
N=1;  % # of Ar
N1=0; % # of Ar+
N2=0; % # of Ar2+
N3=0; % # of Ar3+
N4=0; % # of Ar4+
adk=0;
derN=0;
derN1=0;
derN2=0;
derN3=0;
derN4=0;
Cutoff=0; % Used to Calculate cutoff energy

% Assigned the value for N and omegaADK and cutoff Energy
% Each step is pi/omega1/100*2.42 corresponding to half cycle in fs
% 100fs is the total time range
for i=1:ceil(100*fs/(pi/omega1/100*2.42)) 
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
  
  
  % Based on the cutoff rule Omega_Cutoff = 3.17Up + Ip all the Unit in eV
  % see maxEnergy function 
  % 15.8eV 27.6eV 40.7eV
  % Up = 9.337 38 x 10-5 * I [PW/cm2] ?2 [nm]
  if N(i) >=0.01
      Cutoff0 = maxEnergy('Ne',Inp(2.42*pi/100/omega1*(i-1)-50*fs),lambda);
  else
      Cutoff0 = 0;
  end
  if N1(i)>=0.01
      Cutoff1 = maxEnergy('Ne',Inp(2.42*pi/100/omega1*(i-1)-50*fs),lambda);
  else
      Cutoff1 = 0;
  end 
  if N2(i)>=0.01
      Cutoff2 = maxEnergy('Ne+',Inp(2.42*pi/100/omega1*(i-1)-50*fs),lambda);
  else
      Cutoff2 = 0;
  end 
  if N3(i)>=0.01
      Cutoff3 = maxEnergy('Ne2+',Inp(2.42*pi/100/omega1*(i-1)-50*fs),lambda);
  else
      Cutoff3 = 0;
  end 
  A = [Cutoff0,Cutoff1,Cutoff2,Cutoff3];
  Cutoff(i) = max(A);
  
end

% Block unknown usage
N1(end)=[];
N2(end)=[];
N3(end)=[];
N4(end)=[];
RateN1(end)=[];
RateN2(end)=[];
RateN3(end)=[];
derN1(end)=[];
derN2(end)=[];
derN3(end)=[];
derN4(end)=[];
N(end)=[];
derN(end)=[];

tplot = tslow./100*2.42 ;
% Plot the value of rate and cutoff energy into graph 
figure(1);
plot(tplot,abs(derN));
%legend('Ar','Ar+','Ar2+','The Ar3+ rate*1000')
hold on
plot(tplot,abs(RateN1));
plot(tplot,abs(RateN2));
plot(tplot,1000*abs(RateN3));
legend('Ne','Ne+','Ne2+','The Ne3+ rate*1000')
xlabel('T(fs)');
ylabel('Ionization rate');
yyaxis right;
plot(tslow./100*2.42,Cutoff,'LineWidth',1.5);
ylabel('Cutoff Energy (eV)');
title('NeonPlot')

% Plot the value of fraction and cutoff Energy into graph 
figure(2);
hold on
set(gca,'Yscale','log')
ylim([0.001 1]);
plot(tplot,N,'-c','LineWidth',1);
yyaxis left;
plot(tplot,abs(N1),'-g','LineWidth',1);
plot(tplot,N2,'-r','LineWidth',1);
plot(tplot,N3,'-m','LineWidth',1);
plot(tplot,N4,'-b','LineWidth',1);
grid on
legend('Ne','Ne+','Ne2+','Ne3+','Ne4+');
xlabel('T(fs)');
ylabel('Ionization Fraction');
yyaxis right;
plot(tplot,Cutoff,'-','LineWidth',1.5);
ylabel('Cutoff Energy (eV)');
title('NeonPlot')



