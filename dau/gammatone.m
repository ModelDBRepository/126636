%this is the function to get gammatone filter
%the gain is ~1 when w=w_cf
%
% To use this function, you needn't set the Fs very high(because
% it is bilinear, so just above 2*band is enough and then you may
% need to upsampling the result if you  want to pass the result
% through the synapse ( the synapse dynamics need high sampling rate
%
%Using bilinear transform/Impinvar transform !!!
%function [Numd,Dend] = gammatone(tau,order,cf,Fs) //cf, Fs in Hz
%
function [B,A] = gammatone(tau,order,cf,Fs)

%Z = [];

%K = (1/tau)^order;
w_cf = 2*pi*cf;
a1 = [tau 1-j*tau*w_cf];
a2 = [tau 1+j*tau*w_cf];
a1res = [1];
a2res = [1];

for i = 1:order,
  a1res = conv(a1res,a1);
  a2res = conv(a2res,a2);
end;

num = a1res+a2res;
den = conv(a1res,a2res);
num = [zeros(1,length(den)-length(num)) num];
%[B,A] = bilinear(real(num),real(den),Fs);

[B,A] = impinvar(real(num),real(den),Fs);
