function [ output_args ] = HeliumPlotFunc( lambda,Il0, cycle, it)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


% Clean all the data in previous section

fs=10^(-15)/(2.42*10^(-17));
El0= IE(3.5*10^16);    % Amplitude for the electric field 
                       % Calculated for the laser intensity

c = physconst('LightSpeed');
%Using 258 nm as the wavelength
% omega = 2*pi*c/lamda 
omega1=2*pi*c*10^9*2.42*10^(-17)/lambda;
tau= lambda*10^(-9)*cycle/c/1.76/(2.42*10^(-17));

% 4.3 fs / 1.76
% 4.3 fs is FWHM full width at half maximum 
% get the parameter tau = FWHM/1.76 = 2.44 
E= @(t) El0.*sech(t./tau).*abs(cos(omega1.*t));
Inp = @(t) Il0.*(sech(t./tau).^2); 
%tslow=-50*fs:pi/omega1/100*2.42:50*fs;
tslow=-3*tau*1.76:pi/omega1/it:3*tau*1.76;


% Take each step as pi/omega1/it*2.42 
% need to convert to atomic unit
% Pay attention to the unit conversion
tslow=-50*fs:pi/omega1/it:50*fs;

% adkPeaks=omegaADK(E(tslow),1,.58,1,0); where to find peaks
% define intial value for 
N=1;  % #He
N1=0; % #He+
N2=0; % #He2+
adk=0;
derN=0;
derN1=0;
derN2=0;
Cutoff=0; % Used to Calculate cutoff energy

%ceil rounds the elements of A to the nearest integer
for i=1:ceil(100*fs/(pi/omega1/it)) 
    
  % Calculate out the value of Omega first
  adk(i) = omegaADK(E(pi/it/omega1*(i-1)-50*fs),1,.904,1,0);
  adk1(i)= omegaADK(E(pi/it/omega1*(i-1)-50*fs),2,2,1,0); 
 
  % Calculate the Change of population over time
  % adk multiply the value of each step pi/omega1/it
  % N(i+1)=N(i)*exp(-adk(i)*2.42*pi/omega1/it); %???
  
  % Method to take the derivative to calculate 
  derN(i+1) =adk(i)*N(i);
  derN1(i+1)=adk(i)*N(i)-adk1(i)*N1(i);
  derN2(i+1)=adk1(i)*N1(i);
  
  % Method for using the
  N(i+1) = N(i)-derN(i)*pi/omega1/it;
  N1(i+1)=N1(i)+derN1(i)*pi/omega1/it;
  N2(i+1)=N2(i)+derN2(i)*pi/omega1/it;

  % Final rate should be the omegaADK multiply the value of population
  RateN1(i+1)=adk1(i)*N1(i);
  
  % Based on the cutoff rule Omega_Cutoff = 3.17Up + Ip all the Unit in eV
  % see maxEnergy function 
  % Up = 9.337 38 x 10-5 * I [PW/cm2] ?2 [nm]
  if N(i)>=0.01
      Cutoff0 = maxEnergy('He',Inp(pi/it/omega1*(i-1)-50*fs),lambda);
  else
      Cutoff0 = 0;
  end
  if N1(i)>=0.01
      Cutoff1 = maxEnergy('He+',Inp(pi/it/omega1*(i-1)-50*fs),lambda);
  else
      Cutoff1 = 0;
  end 
  A = [Cutoff0,Cutoff1];
  Cutoff(i) = max(A);
  
end

% Block unknown usage
N1(end)=[];
N2(end)=[];
RateN1(end)=[];
derN1(end)=[];
derN2(end)=[];
N(end)=[];
derN(end)=[];

figure(1);
hold on
yyaxis left ;
plot(tslow./it,abs(derN));
plot(tslow./it,abs(RateN1));
xlabel('T(fs)');
ylabel('Ionization rate');
yyaxis right;
plot(tslow./it,Cutoff);
legend('He','He+','Cutoff Energy');
ylabel('Cutoff Energy (eV)');

figure(2);
hold on
plot(tslow./it,abs(N));
plot(tslow./it,abs(N1));
plot(tslow./it,N2);
xlabel('T(fs)');
ylabel('Ionization Fraction');
yyaxis right;
plot(tslow./it,Cutoff);
legend('He','He+','He2+','Cutoff Energy');
ylabel('Cutoff Energy (eV)');


end

