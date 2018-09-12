
% ==================================================
% Reproduce the data in PRL 103, 143901 (2009)
% Calculation of ionization Rate, ionizatio fraction
% with Cutoff Energy
% ==================================================
clc; clear; % Clean all the data in previous section
fs=10^(-15)/(2.42*10^(-17)); % atomic unit
% Define electric field strength/populations & etc.

Il0 = 3.0376*10^15 ; % Peak intensity Used in Argon
El0 = sqrt(Il0/(1*10^14)) * 0.053376 ;

% 2.42*10^-17 converts the value into atomic unit
% Using laser of wavelength 800nm 
lambda = 800 ;
omega1=2*pi*3*10^8/lambda*10^9*2.42*10^(-17);


% 19fs here used in full pulse 
% tslow calculate the half cycle value
% 10.8 Sech^2 shaped pulses gives that 
% The full width at half-maximum pulse 
% duration is approx 1.76 times the parameter
% E is the electric Field
% Inp is the laser intensity profile 
tau = lambda*10^(-9)/(3*10^8)/(1.76*2.42*10^-17)*7;
E= @(t) El0.*sech(t./tau);
Inp = @(t) Il0.*(sech(t./tau)^2);  
tslow=-50*fs:pi/omega1/100*2.42:50*fs;



% adkPeaks=omegaADK(E(tslow),1,.58,1,0); testing use
% Used to Calculate ionization rate with certain population
N=1;  % # of Ar
N1=0; % # of Ar+
N2=0; % # of Ar2+
N3=0; % # of Ar3+
N4=0; % # of Ar4+
N5=0;
N6=0;
N7=0;

adk=0;
derN=0;
derN1=0;
derN2=0;
derN3=0;
derN4=0;
derN5=0;
derN6=0;
derN7=0;

Cutoff=0; % Used to Calculate cutoff energy


% Assigned the value for N and omegaADK and cutoff Energy
% Each step is pi/omega1/100*2.42 corresponding to half cycle in fs
% 100fs is the total time range
for i=1:ceil(100*fs/(pi/omega1/100*2.42)) 
    %ceil rounds the elements of A to the nearest integer
    
  % Calculate out the value of Omega first
  adk(i)=omegaADK(E(2.42*pi/100/omega1*(i-1)-49.2*fs),1,.58,1,0);
  adk1(i)=omegaADK(E(2.42*pi/100/omega1*(i-1)-49.2*fs),2,1.01,1,0);
  adk2(i)=omegaADK(E(2.42*pi/100/omega1*(i-1)-49.2*fs),3,1.496,1,0);
  adk3(i)=omegaADK(E(2.42*pi/100/omega1*(i-1)-49.2*fs),4,2.1989,1,0);
  
  adk4(i)=omegaADK(E(2.42*pi/100/omega1*(i-1)-49.2*fs),5,2.758,1,0);
  adk5(i)=omegaADK(E(2.42*pi/100/omega1*(i-1)-49.2*fs),3,3.3459,1,0);
  adk6(i)=omegaADK(E(2.42*pi/100/omega1*(i-1)-49.2*fs),4,4.57069853,1,0);
  adk7(i)=omegaADK(E(2.42*pi/100/omega1*(i-1)-49.2*fs),5,5.27426471,1,0);
  
  
  % Calculate the Change of population over time
  % [Might use other code/method to Calculate N]
  N(i+1) = N(i)-derN(i)*2.42*pi/omega1/100;
  % similarly in N1 ... N2 ... N3 ... etc
  %N(i+1)=N(i)*exp(2.42*(-adk(i))*pi/omega1/100);
  
  % Method to take the derivative
  derN(i+1)=adk(i)*N(i);
  derN1(i+1)=adk(i)*N(i)-adk1(i)*N1(i);
  derN2(i+1)=adk1(i)*N1(i)-adk2(i)*N2(i);
  derN3(i+1)=adk2(i)*N2(i)-adk3(i)*N3(i);
  
  derN4(i+1)=adk3(i)*N3(i)-adk4(i)*N4(i);
  derN5(i+1)=adk4(i)*N4(i)-adk5(i)*N5(i);
  derN6(i+1)=adk5(i)*N5(i)-adk6(i)*N6(i);
  derN7(i+1)=adk6(i)*N6(i);
  
  
  % Method for using the 
  N1(i+1)=N1(i)+derN1(i)*2.42*pi/omega1/100;
  N2(i+1)=N2(i)+derN2(i)*2.42*pi/omega1/100;
  N3(i+1)=N3(i)+derN3(i)*2.42*pi/omega1/100;
  N4(i+1)=N4(i)+derN4(i)*2.42*pi/omega1/100;
  N5(i+1)=N5(i)+derN5(i)*2.42*pi/omega1/100;
  N6(i+1)=N6(i)+derN6(i)*2.42*pi/omega1/100;
  N7(i+1)=N7(i)+derN7(i)*2.42*pi/omega1/100;
  
  % Final rate should be the omegaADK multiply the value of population
  RateN1(i+1)=adk1(i)*N1(i);
  RateN2(i+1)=adk2(i)*N2(i);
  RateN3(i+1)=adk3(i)*N3(i);
  RateN4(i+1)=adk4(i)*N4(i);
  RateN5(i+1)=adk5(i)*N5(i);
  RateN6(i+1)=adk6(i)*N6(i);  

  
  % Based on the cutoff rule Omega_Cutoff = 3.17Up + Ip all the Unit in eV
  % see maxEnergy function 
  % Up = 9.337 38 x 10-5 * I [PW/cm2] ?2 [nm]
  if N(i) >=0.01
      Cutoff0 = maxEnergy('Ar',Inp(2.42*pi/100/omega1*(i-1)-49.2*fs),lambda);
      K(i) = Keldysh('Ar',lambda,Inp(2.42*pi/100/omega1*(i-1)-49.2*fs));
  else
      Cutoff0 = 0;
  end
  if N1(i)>=0.01
      Cutoff1 = maxEnergy('Ar',Inp(2.42*pi/100/omega1*(i-1)-49.2*fs),lambda);
      K(i) = Keldysh('Ar',lambda,Inp(2.42*pi/100/omega1*(i-1)-49.2*fs));
  else
      Cutoff1 = 0;
  end 
  if N2(i)>=0.01
      Cutoff2 = maxEnergy('Ar+',Inp(2.42*pi/100/omega1*(i-1)-49.2*fs),lambda);
      K(i) = Keldysh('Ar+',lambda,Inp(2.42*pi/100/omega1*(i-1)-49.2*fs));
  else
      Cutoff2 = 0;
  end 
  if N3(i)>=0.01
      Cutoff3 = maxEnergy('Ar2+',Inp(2.42*pi/100/omega1*(i-1)-49.2*fs),lambda);
      K(i) = Keldysh('Ar2+',lambda,Inp(2.42*pi/100/omega1*(i-1)-49.2*fs));
  else
      Cutoff3 = 0;
  end 
  
  
  if N4(i)>=0.01
      Cutoff4 = maxEnergy('Ar3+',Inp(2.42*pi/100/omega1*(i-1)-49.2*fs),lambda);
      K(i) = Keldysh('Ar3+',lambda,Inp(2.42*pi/100/omega1*(i-1)-49.2*fs));
  else
      Cutoff4 = 0;
  end 
  if N5(i)>=0.01
      Cutoff5 = maxEnergy('Ar4+',Inp(2.42*pi/100/omega1*(i-1)-49.2*fs),lambda);
      K(i) = Keldysh('Ar4+',lambda,Inp(2.42*pi/100/omega1*(i-1)-49.2*fs));
  else
      Cutoff5 = 0;
  end 
  if N6(i)>=0.01
      Cutoff6 = maxEnergy('Ar5+',Inp(2.42*pi/100/omega1*(i-1)-49.2*fs),lambda);
      K(i) = Keldysh('Ar5+',lambda,Inp(2.42*pi/100/omega1*(i-1)-49.2*fs));
  else
      Cutoff6 = 0;
  end 
  if N7(i)>=0.01
      Cutoff7 = maxEnergy('Ar6+',Inp(2.42*pi/100/omega1*(i-1)-49.2*fs),lambda);
      K(i) = Keldysh('Ar6+',lambda,Inp(2.42*pi/100/omega1*(i-1)-49.2*fs));
  else
      Cutoff7 = 0;
  end 
  
  A = [Cutoff0,Cutoff1,Cutoff2,Cutoff3,Cutoff4,Cutoff5,Cutoff6,Cutoff7];
  Cutoff(i) = max(A);
  

  
end

% Block unknown usage
N(end)=[];
N1(end)=[];
N2(end)=[];
N3(end)=[];
N4(end)=[];
N5(end)=[];
N6(end)=[];
N7(end)=[];

RateN1(end)=[];
RateN2(end)=[];
RateN3(end)=[];
RateN4(end)=[];
RateN5(end)=[];
RateN6(end)=[];

derN(end)=[];
derN1(end)=[];
derN2(end)=[];
derN3(end)=[];
derN4(end)=[];
derN5(end)=[];
derN6(end)=[];
derN7(end)=[];




% Plot the value of rate and cutoff energy into graph 
figure(1);
plot(tslow./100*2.42,abs(derN));
%legend('Ar','Ar+','Ar2+','The Ar3+ rate*1000')
hold on
plot(tslow./100*2.42,abs(RateN1));
plot(tslow./100*2.42,abs(RateN2));
plot(tslow./100*2.42,abs(RateN3));
plot(tslow./100*2.42,abs(RateN4));
plot(tslow./100*2.42,abs(RateN5));
plot(tslow./100*2.42,abs(RateN6));
legend('Ar','Ar+','Ar2+','Ar3+','Ar4+','Ar5+','Ar6+')
xlabel('T(fs)');
ylabel('Ionization rate');
yyaxis right;
plot(tslow./100*2.42,Cutoff);
ylabel('Cutoff Energy (eV)');
title('ArgonPlot');


% Plot the value of fraction and cutoff Energy into graph 
figure(2)
fig = figure(2);
left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
hold on
%set(gca,'Yscale','log')
%ylim([0.001 1]);
yyaxis left;
plot(tslow./100*2.42,N,'-c','LineWidth',1);
plot(tslow./100*2.42,abs(N1),'-g','LineWidth',1);
plot(tslow./100*2.42,N2,'-r','LineWidth',1);
plot(tslow./100*2.42,N3,'-m','LineWidth',1);
plot(tslow./100*2.42,N4,'-k','LineWidth',1);
plot(tslow./100*2.42,N5,'LineWidth',1);
plot(tslow./100*2.42,N6,'-b','LineWidth',1);
plot(tslow./100*2.42,N7,'LineWidth',1); 
%plot(tslow./100*2.42,K,'LineWidth',1);
grid on
legend('Ar','Ar+','Ar2+','Ar3+','Ar4+','Ar5+','Ar6+','Ar7+','Keldysh');
xlabel('T(fs)');
ylabel('Ionization Fraction');
yyaxis right;
plot(tslow./100*2.42,Cutoff,'-k','LineWidth',1.5);
ylabel('Cutoff Energy (eV)');
title('ArgonPlot');


