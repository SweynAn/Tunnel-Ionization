
% lowest intensity calculation
% Can you please include ionization by 

% Requirements:
% 400 and 200 nm drivers till depletion of all electrons in He, Ne, Ar?
% Then the same thing for 515, 343, 258 nm drivers?
% Assume 10 cycle drivers at FWHM everywhere. 

% He
I1=LowInt('He+',400,10, 100);
I2=LowInt('He+',200,10, 100);
I3=LowInt('He+',515,10, 100);
I4=LowInt('He+',343,10, 100);
I5=LowInt('He+',258,10, 100);

% Ar
I6=LowInt('Ar5+',400,10, 2000);
I7=LowInt('Ar5+',200,10, 2000);
I8=LowInt('Ar5+',515,10, 2000);
I9=LowInt('Ar5+',343,10, 2000);
I10=LowInt('Ar5+',258,10, 2000);


% Ne
I11=LowInt('Ne2+',400,10, 100);
I12=LowInt('Ne2+',200,10, 100);
I13=LowInt('Ne2+',515,10, 100);
I14=LowInt('Ne2+',343,10, 100);
I15=LowInt('Ne2+',258,10, 100);
