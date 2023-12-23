% code to call the gammatone filter    6/26/02 LHC
% modified by davidson 03/18/03
function filtered = gamma_filt(unfiltered,cf,fs)


% The folowing line specifies Glasberg & Moore style erb's 
erbCF = 24.7 * (4.37 * cf/1000 + 1);

%This converts ERB to a gammatone tau, based on Patterson's
%conversion between Roex and gammatone
tau = 1./(2.*pi * 1.019 * erbCF);

%This calls the fileter in the function gammatone.m
[b1,a1] = gammatone(tau,4,cf,fs);
filtered = filter(b1,a1,unfiltered);

