% Breebaart model code used in Davidson et al, 2009, JASA
clear all
close all

fs = 50e3; %sampling frequency (Hz)
tdres = 1/fs; %time-domain resolution (sec)
dur = 0.3 % duration, seconds
ramplength = 0.01; % ramp duration in seconds (note - duration includes 1/2 of ramp)
f = 500; % target frequency, hz
No = 40; % Noise spectrum level, dB SPL
phi = 0; %tone phase, radians, referenced to sin
To = 20; %target intensity, dB SPL
snr = To - 40;
condition = 2;  % 1 for noso, 2 for nospi;
template_reps = 10; % number of template repetitions
bw = 2900; %noise bandwidth, hz

if bw == 2900  % for Evilsizer et al. data, and Davidson et al 2006 studies
    flf = 275; % low freq of auditory channels used for these stimuli
    fuf = 725;  % high freq of auditory channels used
    fl = 100; %low freq of noise band, hz
    fu = 3000; % high freq of noise band, hz
    sp = 1; % "spacing" of auditory channels (ERBs)
elseif bw == 115; % for Isabelle & Colburn, 1991 stimuli
    flf = 445;
    fuf = 560;
    fl = 445;
    fu = 560;
    sp = .5
elseif bw == 100;  % for Evilsizer et al and Davison et al 2006 studies
    flf = 450;
    fuf = 550;
    fl = 450;
    fu = 550;
    sp= .5
elseif bw == 50; % for Davidson et al. 2009 study
    flf = 475;
    fuf = 525;
    fl = 475;
    fu = 525;
    sp= .2
else
    error('Bad BW')
end
avrms = rms(noise_gen(fs,10,fl,fu)); % compute average rms across a 10 sec noise

N = noise_gen(fs,dur,flf,fuf,No);

S = sin(2*pi*f*[0:tdres:(length(N)/fs)-tdres]+phi);
S = S * sqrt(2) * 20e-6 * 10.0^(To/20); %scaled tone for desired target dB SPL

ramp = hanning(round(dur*fs))';
ramp = [ramp(1:floor(length(ramp)/2)) ones(1,length(S)-length(ramp)) ramp(floor(length(ramp)/2)+1:end) ];
N = (N .* ramp)';
S = (S .* ramp)';

switch condition
    case 1
        R = N+S;
        L = R;
    case 2
        R = N+S;
        L = N-S;
    otherwise
        error('Cannot recognize condition')
end

%------------------------------------------------------------
% Get erb cfs between low and high noise cutoffs
%-----------------------------------------------------------
[finds, fcfs] = getGFBCenterERBs(flf,fuf,sp);
fcfs = erbtofreq(fcfs);

tc = 0; % initialize template counter flag

reps = template_reps+1;

rampts = ramplength * fs;
steadypts = ceil(dur * fs - 2*rampts);
totalpts = steadypts + (rampts*2);

% WAVEFORM ENVELOPE - On/off ramps
step = pi/(rampts-1);
x=[0:step:pi];
offramp = (1+cos(x))./2;
onramp = (1+cos(fliplr(x)))./2;
o=ones(1,(steadypts));
wholeramp = [onramp o offramp]; % Envelope for stimulus (i.e. on/off ramps)

%  tdres = 1/fs;
t=(0:tdres:dur-tdres);

randn('state',sum(100*clock)); % seed the random # generator with the clock
rand('state',sum(100*clock)); % seed the random # generator with the clock

scale_tone = sqrt(2) * 20e-6 * 10.0^(((snr + No)) / 20.0); % in pascals; SPL in dB w.r.t. noise spectral level for 2nd interval

tone_ramp = wholeramp;

tone=sin(2*pi*f*t);

rmsnoise = No + (10 .* log10(bw));
scale_noise = 20e-6*10^((rmsnoise)/20);  % into Pascals

lofreq = fl;
bw = fu-fl;
hifreq = fu;

for rep = 1:reps
    totalpts = length(R);
    stim1 = randn(1,totalpts);
    stim2 = randn(1,totalpts);
    fres = 1./(totalpts * tdres);
    
    nptslofreq = floor(lofreq/fres);
    nptshifreq = ceil(hifreq/fres);
    % Now filter the noise
    stimspectrum = fft(stim1);
    stimspectrum(1:nptslofreq) = 0.;  % zero out all bands that apply (Mickey's bug-fix in here)
    stimspectrum(nptshifreq:(totalpts - nptshifreq+2)) = 0.;
    stimspectrum((totalpts - nptslofreq+2):totalpts) = 0.;
    % So, values in the fft of stim1 (both real & imag parts) have been zeroed
    stim1a = real(ifft(stimspectrum)); % inverse FFT to create noise waveform
    
    stimspectrum = fft(stim2);
    stimspectrum(1:nptslofreq) = 0.;  % zero out all bands that apply (Mickey's bug-fix in here)
    stimspectrum(nptshifreq:(totalpts - nptshifreq+2)) = 0.;
    stimspectrum((totalpts - nptslofreq+2):totalpts) = 0.;
    % So, values in the fft of stim1 (both real & imag parts) have been zeroed
    stim2a = real(ifft(stimspectrum)); % inverse FFT to create noise waveform
    
    if tc == 0;
        stim1 = ((stim1a/rms(stim1a) .* wholeramp * scale_noise)+ (scale_tone * tone .* tone_ramp)) * 5e4;
        stim2 = ((stim1a/rms(stim1a) .* wholeramp * scale_noise)+ (-scale_tone * tone .* tone_ramp)) * 5e4;
        stim3 = ((stim2a/rms(stim2a)) .* wholeramp * scale_noise)* 5e4;
    else
        stim1 = R(:,rep-template_reps)' * 5e4;
        stim2 = L(:,rep-template_reps)' * 5e4;
        stim3 = N(:,rep-template_reps)' * 5e4;
    end
    clear pstim1 pstim2 pstim1f pstim2f fcfs;
    
    [fcfs, pstim1]= breb_periph(fs,fl,fu,stim1,sp);
    [fcfs, pstim2]= breb_periph(fs,fl,fu,stim2,sp);
    [fcfs, pstim3]= breb_periph(fs,fl,fu,stim3,sp);
    
    c = 0.030; % constant from Breebaart et al, 2001
    
    num_f = [ (1-exp(-1/(fs*c))) 0];
    den_f = [ 1 -exp(-1/(fs*c))];
    pst_scl = 1-exp(-1/(fs*c));
    
    pstimS = (pstim2-pstim1).^2;
    pstimN = (pstim3-pstim3).^2;
    
    % filter using exponential
    pstimS = filter(num_f,den_f,pstimS)+ flipud(filter(num_f,den_f,flipud(pstimS)))- (pst_scl * pstimS);
    pstimN = filter(num_f,den_f,pstimN)+ flipud(filter(num_f,den_f,flipud(pstimN)))- (pst_scl * pstimN);
    
    a = 0.1; b = 0.00002; % constant from Breebaart et al, 2001
    pstimS = a*log10(b*pstimS+1); %apply nonlinear transformation
    pstimN = a*log10(b*pstimN+1);
    
    if  ((tc==0) & (template_reps > 0));
        template_s(rep,:,:) = pstimS';
        template(rep,:,:) = pstimN';
        disp('Template Refresh')
        if (rep==template_reps);
            tc = 1;
        end
    else
        %apply noise if desired, to come up with E prime
        E_prime2(1,:,:) = pstimS'; % + noise_sd*randn(1,length(pstimS)));
        E_prime2_noise(1,:,:) = pstimN';% + noise_sd*randn(1,length(pstimN)));
    end
    
    if rep>template_reps;
        Ep2 = squeeze(E_prime2);
        Ep2n = squeeze(E_prime2_noise);
        
        template_c = squeeze(mean(template));
        template_cat = template;
        
        var_n = 1 + squeeze(mean((template_cat.^2)))-template_c.^2;
        
        template_s_mu = squeeze(mean(template_s,1));
        mu = template_s_mu - template_c;
        
        %come up with weights
        int_n_sigma = sqrt(sum(sum(squeeze((mu.^2)./(var_n.^2)))));
        
        %compute weighted distance from template and sum across time and freq
        Us(rep-template_reps) = sum(sum(((mu./var_n).*(Ep2-template_c))));
        Un(rep-template_reps) = sum(sum(((mu./var_n).*(Ep2n-template_c))));
    end
end

R2 = Us % decision variables (see below)
N2 = Un

% Note: dp = inf, unless internal noise is added!!! (not included in Davidson et al 2009)
dp = (mean(R2)-mean(N2)) / sqrt(mean([var(N2) var(R2)])) % dprime

%In Davidson et al, 2009, model outputs were correlated to data, as follows:
%         py = [data.pd{i}(:,s) ; data.pf{i}(:,s)]; % Probability of "yes" is the combined data for prob detection and prob false alarms 
%         zy = zscore(py);
%         corr_py{idata_set}(i,s,:) = min(min(corrcoef(zy,[R2 N2]))); % compute correlation between zy and concatenated model outputs for T+N and N
%         corr_pd{idata_set}(i,s,:) = min(min(corrcoef(zscore(data.pd{i}(:,s)),R2))); % correlation  between prob detection and model responses for T+N stimuli
%         corr_pf{idata_set}(i,s,:) = min(min(corrcoef(zscore(data.pf{i}(:,s)),N2))); % correlation between prob FA and Model output for N stimuli
%         
%         dv{idata_set}(i,s,:) = [R2 N2]; % decision variables for model output

        
          
   