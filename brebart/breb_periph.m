% Peripheral Model for Breebaart model
%input in pascals
%output in model units

function  [fcfs,out ] = breb_periph(fs,fl,fu,stim,fspace);
tdres = 1/fs;

if nargin <5;
[finds, fcfs] = getGFBCenterERBs(fl,fu,.5);
else 
  if (fl ==100) & (fu==3000)  
    [finds, fcfs] = getGFBCenterERBs(275,725,fspace);  
  else
    [finds, fcfs] = getGFBCenterERBs(fl,fu,fspace); 
  end  
end
fcfs = erbtofreq(fcfs);
       %%%%% outer and middle ear transfer fcn
                q = 2-cos(2*pi*4000/fs)-sqrt(((cos(2*pi*4000/fs))-2)^2-1);
                r = 2-cos(2*pi*1000/fs)-sqrt(((cos(2*pi*1000/fs))-2)^2-1);
                num_mi = [(1-q)*r -(1-q)*r 0];
                den_mi = [1 -(q+r) q*r];
                
                stim2 = filter(num_mi,den_mi,stim);
                
                %%%%% Basilar Membrane filtering
                v = 3;

                n = 0;
                for fc = fcfs
                    n = n+1;
                    stim3 = stim2.*exp(-2*pi*i*fc*([1:length(stim2)]*tdres));
                    b = (24.7*(0.00437*fc+1)*factorial(v-1)^2)/(pi*factorial(2*v-2)*2^-(2*v-2));
                    num_gt = [(1-exp(-2*pi*b*tdres)) 0];
                    den_gt = [1 -exp(-2*pi*b*tdres)];
                    for x = 1:3;
                       stim3 = filter(num_gt,den_gt,stim3);
                    end
                    stim4(:,n) = 2*real(stim3.*exp(2*pi*i*fc*([1:length(stim3)]*tdres)))';
                end
                %internal noise
                stim5 = stim4+ (randn(size(stim4))*60e-6);
       
                %half-wave rect
                stim6 = stim5.*(stim5>0);
         
                %loss of synch
                u = 2-cos(2*pi*2000/fs)-sqrt(((cos(2*pi*2000/fs)-2)^2)-1);
                num_pl = [1-u 0];
                den_pl = [1 -u];
                stim6b = stim6;
         
                for x = 1:5
                    stim6b = filter(num_pl,den_pl,stim6b);
                end
  
                %nonlinear adaptation loops
            for n = 1:length(fcfs)
                stim8(:,n) = fadapt(stim6b(:,n),fs);
            end
       out = stim8;
   
                