
% ==================================================
% Reproduce the data in PRL 103, 143901 (2009)
% Calculation of ionization Rate, ionizatio fraction
% with Cutoff Energy
% ==================================================
clc; clear; % Clean all the data in previous section
fs=10^(-15)/(2.42*10^(-17)); % atomic unit
% Define electric field strength/populations & etc.

Il0 = 2.1120e+18 ; % Peak intensity Used in Argon
El0 = sqrt(Il0/(1*10^14)) * 0.053376 ;

% 2.42*10^-17 converts the value into atomic unit
% Using laser of wavelength 800nm 
lambda = 500 ;
omega1=2*pi*3*10^8/lambda*10^9*2.42*10^(-17);
it = 200000;
cycle = 10;

% 19fs here used in full pulse 
% tslow calculate the half cycle value
% 10.8 Sech^2 shaped pulses gives that 
% The full width at half-maximum pulse 
% duration is approx 1.76 times the parameter
% E is the electric Field
% Inp is the laser intensity profile 
tau = lambda*10^(-9)/(3*10^8)/(1.76*2.42*10^-17)*cycle;
E= @(t) El0.*sech(t./tau).*abs(cos(omega1.*t));
Inp = @(t) Il0.*(sech(t./tau)^2);  
tslow=-50*fs:pi/omega1/it*2.42:50*fs;


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

Cutoff=0; % Used to Calculate cutoff energy


% Assigned the value for N and omegaADK and cutoff Energy
% Each step is pi/omega1/it*2.42 corresponding to half cycle in fs
% itfs is the total time range
for i=1:ceil(100*fs/(pi/omega1/it*2.42)) 
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
  derN1(i+1)=adk(i)*N(i)  - adk1(i)*N1(i);
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
  
  % Final rate should be the omegaADK multiply the value of population
  RateN1(i+1)=adk1(i)*N1(i);
  RateN2(i+1)=adk2(i)*N2(i);
  RateN3(i+1)=adk3(i)*N3(i);
  RateN4(i+1)=adk4(i)*N4(i);
  RateN5(i+1)=adk5(i)*N5(i);
  RateN6(i+1)=adk6(i)*N6(i);  
  % Update to the rate
  RateN7(i+1)=adk7(i)*N7(i);
  RateN8(i+1)=adk8(i)*N8(i);
  RateN9(i+1)=adk9(i)*N9(i);  
  RateN10(i+1)=adk10(i)*N10(i);
  RateN11(i+1)=adk11(i)*N11(i);
  RateN12(i+1)=adk12(i)*N12(i);  
  RateN13(i+1)=adk13(i)*N13(i);
  RateN14(i+1)=adk14(i)*N14(i);
  RateN15(i+1)=adk15(i)*N15(i);
  RateN16(i+1)=adk16(i)*N16(i);
  % Update to higher ions
  
  % Based on the cutoff rule Omega_Cutoff = 3.17Up + Ip all the Unit in eV
  % see maxEnergy function 
  % Up = 9.337 38 x 10-5 * I [PW/cm2] ?2 [nm]
  if N(i) >=0.01
      Cutoff0 = maxEnergy('Ar',Inp(2.42*pi/it/omega1*(i-1)-50*fs),lambda);
      K(i) = Keldysh('Ar',lambda,Inp(2.42*pi/it/omega1*(i-1)-50*fs));
  else
      Cutoff0 = 0;
  end
  
  if N1(i)>=0.01
      Cutoff1 = maxEnergy('Ar',Inp(2.42*pi/it/omega1*(i-1)-50*fs),lambda);
      K(i) =    Keldysh('Ar',lambda,Inp(2.42*pi/it/omega1*(i-1)-50*fs));
  else
      Cutoff1 = 0;
  end 
  if N2(i)>=0.01
      Cutoff2 = maxEnergy('Ar+',Inp(2.42*pi/it/omega1*(i-1)-50*fs),lambda);
      K(i) = Keldysh('Ar+',lambda,Inp(2.42*pi/it/omega1*(i-1)-50*fs));
  else
      Cutoff2 = 0;
  end 
  if N3(i)>=0.01
      Cutoff3 = maxEnergy('Ar2+',Inp(2.42*pi/it/omega1*(i-1)-50*fs),lambda);
      K(i) = Keldysh('Ar2+',lambda,Inp(2.42*pi/it/omega1*(i-1)-50*fs));
  else
      Cutoff3 = 0;
  end 
  
  
  if N4(i)>=0.01
      Cutoff4 = maxEnergy('Ar3+',Inp(2.42*pi/it/omega1*(i-1)-50*fs),lambda);
      K(i) = Keldysh('Ar3+',lambda,Inp(2.42*pi/it/omega1*(i-1)-50*fs));
  else
      Cutoff4 = 0;
  end 
  if N5(i)>=0.01
      Cutoff5 = maxEnergy('Ar4+',Inp(2.42*pi/it/omega1*(i-1)-50*fs),lambda);
      K(i) = Keldysh('Ar4+',lambda,Inp(2.42*pi/it/omega1*(i-1)-50*fs));
  else
      Cutoff5 = 0;
  end 
  if N6(i)>=0.01
      Cutoff6 = maxEnergy('Ar5+',Inp(2.42*pi/it/omega1*(i-1)-50*fs),lambda);
      K(i) = Keldysh('Ar5+',lambda,Inp(2.42*pi/it/omega1*(i-1)-50*fs));
  else
      Cutoff6 = 0;
  end 
  if N7(i)>=0.01
      Cutoff7 = maxEnergy('Ar6+',Inp(2.42*pi/it/omega1*(i-1)-50*fs),lambda);
      K(i) = Keldysh('Ar6+',lambda,Inp(2.42*pi/it/omega1*(i-1)-50*fs));
  else
      Cutoff7 = 0;
  end
  
  
  % update to higher ionization
  if N8(i)>=0.01
      Cutoff8 = maxEnergy('Ar7+',Inp(2.42*pi/it/omega1*(i-1)-50*fs),lambda);
      K(i) = Keldysh('Ar7+',lambda,Inp(2.42*pi/it/omega1*(i-1)-50*fs));
  else
      Cutoff8 = 0;
  end 
  if N9(i)>=0.01
      Cutoff9 = maxEnergy('Ar8+',Inp(2.42*pi/it/omega1*(i-1)-50*fs),lambda);
      K(i) = Keldysh('Ar8+',lambda,Inp(2.42*pi/it/omega1*(i-1)-50*fs));
  else
      Cutoff9 = 0;
  end 
 
  if N10(i)>=0.01
      Cutoff10 = maxEnergy('Ar9+',Inp(2.42*pi/it/omega1*(i-1)-50*fs),lambda);
      K(i) = Keldysh('Ar9+',lambda,Inp(2.42*pi/it/omega1*(i-1)-50*fs));
  else
      Cutoff10 = 0;
  end 
  
  if N11(i)>=0.01
      Cutoff11 = maxEnergy('Ar10+',Inp(2.42*pi/it/omega1*(i-1)-50*fs),lambda);
      K(i) = Keldysh('Ar10+',lambda,Inp(2.42*pi/it/omega1*(i-1)-50*fs));
  else
      Cutoff11 = 0;
  end 
  
  if N12(i)>=0.01
      Cutoff12 = maxEnergy('Ar11+',Inp(2.42*pi/it/omega1*(i-1)-50*fs),lambda);
      K(i) = Keldysh('Ar11+',lambda,Inp(2.42*pi/it/omega1*(i-1)-50*fs));
  else
      Cutoff12 = 0;
  end 
  
  if N13(i)>=0.01
      Cutoff13 = maxEnergy('Ar12+',Inp(2.42*pi/it/omega1*(i-1)-50*fs),lambda);
      K(i) = Keldysh('Ar12+',lambda,Inp(2.42*pi/it/omega1*(i-1)-50*fs));
  else
      Cutoff13 = 0;
  end 

  
  if N14(i)>=0.01
      Cutoff14 = maxEnergy('Ar13+',Inp(2.42*pi/it/omega1*(i-1)-50*fs),lambda);
      K(i) = Keldysh('Ar13+',lambda,Inp(2.42*pi/it/omega1*(i-1)-50*fs));
  else
      Cutoff14 = 0;
  end 
  
  if N15(i)>=0.01
      Cutoff15 = maxEnergy('Ar14+',Inp(2.42*pi/it/omega1*(i-1)-50*fs),lambda);
      K(i) = Keldysh('Ar14+',lambda,Inp(2.42*pi/it/omega1*(i-1)-50*fs));
  else
      Cutoff15 = 0;
  end   
  
  %=====
  if N16(i)>=0.01
      Cutoff16 = maxEnergy('Ar15+',Inp(2.42*pi/it/omega1*(i-1)-50*fs),lambda);
      K(i) = Keldysh('Ar15+',lambda,Inp(2.42*pi/it/omega1*(i-1)-50*fs));
  else
      Cutoff16 = 0;
  end   
  
  if N17(i)>=0.01
      Cutoff17 = maxEnergy('Ar16+',Inp(2.42*pi/it/omega1*(i-1)-50*fs),lambda);
      K(i) = Keldysh('Ar16+',lambda,Inp(2.42*pi/it/omega1*(i-1)-50*fs));
  else
      Cutoff17 = 0;
  end   
  
  
  A = [Cutoff0,Cutoff1,Cutoff2,Cutoff3,Cutoff4,Cutoff5,Cutoff6,Cutoff7, Cutoff8,...
      Cutoff9, Cutoff10, Cutoff11, Cutoff12, Cutoff13, Cutoff14, Cutoff15, Cutoff16, Cutoff17];
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
N8(end)=[];
N9(end)=[];
N10(end)=[];
N11(end)=[];
N12(end)=[];
N13(end)=[];
N14(end)=[];
N15(end)=[];
N16(end)=[];
N17(end)=[];

RateN1(end)=[];
RateN2(end)=[];
RateN3(end)=[];
RateN4(end)=[];
RateN5(end)=[];
RateN6(end)=[];
RateN7(end)=[];
RateN8(end)=[];
RateN9(end)=[];
RateN10(end)=[];
RateN11(end)=[];
RateN12(end)=[];
RateN13(end)=[];
RateN14(end)=[];
RateN15(end)=[];
RateN16(end)=[];


derN(end)=[];
derN1(end)=[];
derN2(end)=[];
derN3(end)=[];
derN4(end)=[];
derN5(end)=[];
derN6(end)=[];
derN7(end)=[];
derN8(end)=[];
derN9(end)=[];
derN10(end)=[];
derN11(end)=[];
derN12(end)=[];
derN13(end)=[];
derN14(end)=[];
derN15(end)=[];
derN16(end)=[];
derN17(end)=[];

% Plot the value of rate and cutoff energy into graph 
figure(1);
plot(tslow./it*2.42,abs(derN));
%legend('Ar','Ar+','Ar2+','The Ar3+ rate*it0')
hold on
plot(tslow./it*2.42,abs(RateN1));
plot(tslow./it*2.42,abs(RateN2));
plot(tslow./it*2.42,abs(RateN3));
plot(tslow./it*2.42,abs(RateN4));
plot(tslow./it*2.42,abs(RateN5));
plot(tslow./it*2.42,abs(RateN6));
plot(tslow./it*2.42,abs(RateN7));
plot(tslow./it*2.42,abs(RateN8));
plot(tslow./it*2.42,abs(RateN9));
plot(tslow./it*2.42,abs(RateN10));
plot(tslow./it*2.42,abs(RateN11));
plot(tslow./it*2.42,abs(RateN12));
legend('Ar','Ar+','Ar2+','Ar3+','Ar4+','Ar5+','Ar6+')
xlabel('T(fs)');
ylabel('Ionization rate');
yyaxis right;
plot(tslow./it*2.42,Cutoff);
ylabel('Cutoff Energy (eV)');
title('ArgonPlot');


% Plot the value of fraction and cutoff Energy into graph 
figure(2)
fig = figure(2);
left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
hold on
set(gca,'Yscale','log')
ylim([0.001 1]);
yyaxis left;
plot(tslow./it*2.42,N,'-c','LineWidth',1);
plot(tslow./it*2.42,abs(N1),'-g','LineWidth',1);
plot(tslow./it*2.42,N2,'-r','LineWidth',1);
plot(tslow./it*2.42,N3,'-m','LineWidth',1);
plot(tslow./it*2.42,N4,'-k','LineWidth',1);
plot(tslow./it*2.42,N5,'LineWidth',1);
plot(tslow./it*2.42,N6,'-b','LineWidth',1);
plot(tslow./it*2.42,N7,'LineWidth',1); 
plot(tslow./it*2.42,N8,'LineWidth',1);
plot(tslow./it*2.42,N9,'-b','LineWidth',1);
plot(tslow./it*2.42,N10,'LineWidth',1); 
plot(tslow./it*2.42,abs(N11),'-g','LineWidth',1);
plot(tslow./it*2.42,N12,'-r','LineWidth',1);
plot(tslow./it*2.42,N13,'-m','LineWidth',1);
plot(tslow./it*2.42,N14,'-k','LineWidth',1);
plot(tslow./it*2.42,N15,'LineWidth',1);
plot(tslow./it*2.42,N16,'-b','LineWidth',1);
plot(tslow./it*2.42,N17,'LineWidth',1); 

%plot(tslow./it*2.42,K,'LineWidth',1);
grid on
legend('Ar','Ar+','Ar2+','Ar3+','Ar4+','Ar5+','Ar6+','Ar7+','Keldysh');
xlabel('T(fs)');
ylabel('Ionization Fraction');
yyaxis right;
plot(tslow./it*2.42,Cutoff,'-k','LineWidth',1.5);
ylabel('Cutoff Energy (eV)');
title('ArgonPlot');


