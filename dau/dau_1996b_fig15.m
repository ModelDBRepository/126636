
close all
clear all
 
%This is a 300-line m file designed to unify the code from Dau's 1996 model and the PEMO/AFC package for ultra-fast implementation.


%%%%%%%%%%%%%%%%%%%%  Stimulus Parameters %%%%%%%%%%%%%%%%%%%%%
fs = 50000;  %sampling rate
tdres = 1/fs; % time-doimain resolution

bw = 4980;   %bandwidth of stimulus noise 
lofreq = 20; % low frequency bound of stimulus noise
No = 40;     % spectrum level of noise
 tmp = 0.3410; %average RMS of this noise over a few reps. This should be recomputed when changing noise parameters
    
snr  =  45;    % signal-to-noise ratio (dB SPL above spectrum level of noise) NOT TONE LEVEL
centerfreq = 1000;  % Tone frequency
dur = 0.300; % stimulus duration in seconds 
ramplength = 0.01; % onset/offset ramp duration



%%%%%%%%%%%%%%%%%% Model parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
template_snr = 45; % dB above noise spectrum level for template 
template_reps = 10; % number of reps used for template
n_track = 1; %number of tracks to run
step_size = [4 2 1 0.5]; %step sizes
step_revs = [2 2 2 6]; % Number of reversals to complete at each step size. 
                       % Thresh is the median of the reversals at the last step size



%%%%%%%%%%%%%%%%%%%%%%%%%% Stimulus Construction

hifreq = bw+lofreq;
rampts = ramplength * fs;
steadypts = ceil(dur * fs - 2*rampts);
totalpts = steadypts + (rampts*2);

step = pi/(rampts-1);
x=[0:step:pi];
offramp = (1+cos(x))./2;
onramp = (1+cos(fliplr(x)))./2;
o=ones(1,(steadypts));
wholeramp = [onramp o offramp]; % Envelope for stimulus (i.e. on/off ramps)   

onenvdur = 0.13;
tonedur = 0.005;
offenvdur = dur-onenvdur-tonedur;
tone_ramp = [ zeros(1,onenvdur*fs) hanning(floor(tonedur*fs))' 0 zeros(1,floor(offenvdur*fs)) ]; % Envelope for stimulus (i.e. on/off ramps)


t=(0:tdres:dur-tdres); 

randn('state',sum(100*clock)); % seed the random # generator with the clock
rand('state',sum(100*clock)); % seed the random # generator with the clock
scale_tone_templ = sqrt(2) * 20e-6 * 10.0^(((template_snr + No)) / 20.0);
tone=sin(2*pi*centerfreq*t);
tone180 = sin(2*pi*centerfreq*t + pi);


rmsnoise = No + (10 .* log10(bw));   
scale_noise = 20e-6*10^((rmsnoise)/20);  % into Pascals




%%%%%%%%%%%%%%%%%%%%% Low-pass filtering coefficients
[b8,a8] = IRIfolp(8,fs); % to achieve a 20 msec time constant (enevlope extraction)
[b1000,a1000] = IRIfolp(1000,fs); % remove fine structure



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Model Code
for i_track = 1:n_track;
    
    template_n_mat = []; %Initalize noise-alone template
    template_tin_mat = []; %and tone-plus-noise template 
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%CREATE TEMPLATE    
    disp(['Computing Template: Track ' num2str(i_track) ])
    for i = 1:template_reps
        % disp(['Template REP ' num2str(i)])
        stim1 = randn(1,totalpts);  
        stim2 = randn(1,totalpts);
        fres = 1./(totalpts * tdres);
        
        % Now filter the noise 
        stimspectrum = fft(stim1);
        nptslofreq = floor(lofreq/fres);
        nptshifreq = ceil(hifreq/fres);
        stimspectrum(1:nptslofreq) = 0.;  % zero out all bands that apply 
        stimspectrum(nptshifreq:(totalpts - nptshifreq+2)) = 0.;
        stimspectrum((totalpts - nptslofreq+2):totalpts) = 0.;
        % So, values in the fft of stim1 (both real & imag parts) have been zeroed
        stim1 = real(ifft(stimspectrum)); % inverse FFT to create noise waveform
        
        stimspectrum = fft(stim2);
        nptslofreq = floor(lofreq/fres);
        nptshifreq = ceil(hifreq/fres);
        stimspectrum(1:nptslofreq) = 0.;  % zero out all bands that apply (Mickey's bug-fix in here)
        stimspectrum(nptshifreq:(totalpts - nptshifreq+2)) = 0.;
        stimspectrum((totalpts - nptslofreq+2):totalpts) = 0.;
        % So, values in the fft of stim1 (both real & imag parts) have been zeroed
        stim2 = real(ifft(stimspectrum)); % inverse FFT to create noise waveform
        
        template_n = stim1/tmp .* wholeramp * scale_noise ; %normalize by rms and scale up to desired rms   
        template_tin = (stim2/tmp .* wholeramp * scale_noise)+ (scale_tone_templ * (tone .* tone_ramp));%/rms(tone.*tone_ramp));  % Now add tone to for "hit" trials  ( but not for FA trials)
        
        template_fscaled_n = gamma_filt(template_n,centerfreq,fs); %gammatone filter
        template_fscaled_tin = gamma_filt(template_tin,centerfreq,fs);
        
        template_fscaled_n = template_fscaled_n .* (template_fscaled_n>0);% half-wave rectify
        template_fscaled_tin = template_fscaled_tin .* (template_fscaled_tin>0);
        
        template_fscaled2_n = filter(b1000,a1000,template_fscaled_n); %low pass filter at 1000 Hz
        template_fscaled2_tin = filter(b1000,a1000,template_fscaled_tin);
        
        
        template_adap_n = nlalmex(template_fscaled2_n',fs)'; %adaptation loops
        template_adap_tin = nlalmex(template_fscaled2_tin',fs)';
        
        template_fadap_n = filter(b8,a8,template_adap_n); %lowpass filter again at 8 Hz to get envelope
        template_fadap_tin = filter(b8,a8,template_adap_tin);
        
        template_n_mat = [template_n_mat; template_fadap_n];    %Store noise alone template
        template_tin_mat = [template_tin_mat; template_fadap_tin]; %Store T+N template
    end
    
    
    template_n = mean(template_n_mat);   %compute mean of templates
    template_tin = mean(template_tin_mat);
    
    template_ndiff = (template_tin - template_n) ; %average difference template
    template_ndiff = template_ndiff/rms(template_ndiff); % compute the normalized difference template
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%PROCESS STIMULI (three interval forced choice)
    
    
    %tracking stuff;
    
    cur_revs = 0;
    direct = -1;
    prev_correct = 0;
    snr_revs =[];
    snr_trials = [];
    direct_trials = [];
    i_step_rev = 1;
    just_change = 0;
    
    
    while cur_revs ~= step_revs(i_step_rev);
        
         disp(['Tracking: Tone level ' num2str(snr+No) 'dB SPL'])
        
        %recompute scaletone for current SNR
        
        scale_tone = sqrt(2) * 20e-6 * 10.0^(((snr + No)) / 20.0); % 
        
        stim1 = randn(1,totalpts);  
        stim2 = randn(1,totalpts);
        stim3 = randn(1,totalpts);
        fres = 1./(totalpts * tdres);
        
        % Now filter the noise 
        stimspectrum = fft(stim1);
        nptslofreq = floor(lofreq/fres);
        nptshifreq = ceil(hifreq/fres);
        stimspectrum(1:nptslofreq) = 0.;  % zero out all bands that apply (Mickey's bug-fix in here)
        stimspectrum(nptshifreq:(totalpts - nptshifreq+2)) = 0.;
        stimspectrum((totalpts - nptslofreq+2):totalpts) = 0.;
        % So, values in the fft of stim1 (both real & imag parts) have been zeroed
        stim1 = real(ifft(stimspectrum)); % inverse FFT to create noise waveform
        stim1 = stim1/tmp * scale_noise .* wholeramp; %normalize by rms and scale up to desired rms   
        
        stimspectrum = fft(stim2);
        nptslofreq = floor(lofreq/fres);
        nptshifreq = ceil(hifreq/fres);
        stimspectrum(1:nptslofreq) = 0.;  % zero out all bands that apply (Mickey's bug-fix in here)
        stimspectrum(nptshifreq:(totalpts - nptshifreq+2)) = 0.;
        stimspectrum((totalpts - nptslofreq+2):totalpts) = 0.;
        % So, values in the fft of stim1 (both real & imag parts) have been zeroed
        stim2 = real(ifft(stimspectrum)); % inverse FFT to create noise waveform
        stim2 = stim2/tmp * scale_noise .* wholeramp; %normalize by rms and scale up to desired rms   
        
        
        %ADD Tone to stim 2
        stim2 = (stim2 )+   (scale_tone * (tone .* tone_ramp));%/rms(tone.*tone_ramp));  % Now add tone to for "hit" trials  ( but not for FA trials)
        
        stimspectrum = fft(stim3);
        nptslofreq = floor(lofreq/fres);
        nptshifreq = ceil(hifreq/fres);
        stimspectrum(1:nptslofreq) = 0.;  % zero out all bands that apply (Mickey's bug-fix in here)
        stimspectrum(nptshifreq:(totalpts - nptshifreq+2)) = 0.;
        stimspectrum((totalpts - nptslofreq+2):totalpts) = 0.;
        % So, values in the fft of stim1 (both real & imag parts) have been zeroed
        stim3 = real(ifft(stimspectrum)); % inverse FFT to create noise waveform
        stim3 = stim3/tmp * scale_noise .* wholeramp; %normalize by rms and scale up to desired rms   
        
        
        
        %filter scaled stimuli
        noise_fscaled = gamma_filt(stim1,centerfreq,fs);  %noise-alone #1
        tin_fscaled = gamma_filt(stim2,centerfreq,fs);    %tone plus noise
        noise3_fscaled = gamma_filt(stim3,centerfreq,fs); %noise-alone #2
        
        %half-wave rectify
        noise_fscaled = noise_fscaled .* (noise_fscaled>0);
        tin_fscaled = tin_fscaled .* (tin_fscaled>0);
        noise3_fscaled = noise3_fscaled .* (noise3_fscaled>0);
        
        %low-pass filter
        noise_fscaled2 = filter(b1000,a1000,noise_fscaled);
        tin_fscaled2 = filter(b1000,a1000,tin_fscaled);
        noise3_fscaled2 = filter(b1000,a1000,noise3_fscaled);
        
        %non-linear adaptation loops
        noise_adap = nlalmex(noise_fscaled2',fs)';
        tin_adap = nlalmex(tin_fscaled2',fs)';
        noise3_adap = nlalmex(noise3_fscaled2',fs)';
        
        %low-pass filter Tau = 20msec
        noise_fadap_intn = filter(b8,a8,noise_adap); 
        tin_fadap_intn = filter(b8,a8,tin_adap);
        noise3_fadap_intn = filter(b8,a8,noise3_adap);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Optimal Detector
        %The optimal detector first computes the difference between the current
        %internal representation on each interval and the average noise-alone template (the noise-alone template was computed above), and then computes
        %the salar product of this difference and the normalized difference template (the normalized-difference template was also computed above) 
        
        % The scalar product is the result of the derivation in the appendix of the 1996a paper
        % The normalization effectively gets rid of the second term which is related to energy
        corrDTU_noise_dtmp = mean((noise_fadap_intn-template_n).* template_ndiff);
        corrDTU_noise3_dtmp = mean((noise3_fadap_intn-template_n).* template_ndiff);
        corrDTU_tin_dtmp =  mean((tin_fadap_intn-template_n).*template_ndiff);
        
        %This step was taken from Torstens code and is a clever way of incoporating internal noise and computing a probability of getting a correct 
        %answer on each trial based on the difference between the tone-plus-noise scalar product and the larger of the noise-alone scalar products
        % the constants adjust the psychometric based on a 2,3,4 interval task etc.  I've included the original code as well. It's in a file called pemo_mdecide.m
        %     
        
        prob = 1 - (erfc(((((corrDTU_tin_dtmp -  max([corrDTU_noise3_dtmp corrDTU_noise_dtmp])) / sqrt(3)) * 0.765) - 0.423) * 0.7071068) / 2); %3-int
        % 
        
        %make the decision
        if prob > rand
            correct = 1;
        else
            correct = 0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tracking Logic %%%%%%%%%%%%%%%%%%%%%%%%%
        snr_trials = [snr_trials snr];
        direct_trials = [direct_trials direct];
        if ((correct == 1) & (prev_correct == 1)); % 2 correct in a row
            
            if direct == 1; %if going up
                direct = -1;%then reversal
                cur_revs = cur_revs + 1; % update reversal counter
                snr_revs = [snr_revs snr]; % log snr;
                snr = snr - step_size(i_step_rev); % lower snr by stepsize  
            elseif direct == -1; %then going down (no reversal)
                
                if snr_trials(end-1) == snr_trials(end); % if the current snr == the previous snr then it's time to lower it; 
                    snr = snr - step_size(i_step_rev); % lower snr by stepsize
                end
            end
            
            prev_correct = 1; 
            
            
        elseif ((correct == 1) & (prev_correct == 0)); % 1 correct 
            prev_correct = 1;  %indicate that this trial was correct;
            
        elseif (correct == 0);  % 1 wrong
            if direct == -1; % if going down;
                direct = 1; %then reversal;
                cur_revs = cur_revs + 1; % update reversal counter
                snr_revs = [snr_revs snr]; % log snr;
                snr = snr + step_size(i_step_rev); %  increase snr by stepsize
            elseif direct == 1; %if going up (no reversal)
                snr_revs = [snr_revs snr]; % log snr;
                snr = snr + step_size(i_step_rev); %  increase snr by stepsize   
            end    
            prev_correct = 0; %indicate that this trial was wrong!
        end 
        
        
        if cur_revs == step_revs(i_step_rev) 
            if i_step_rev == length(step_revs);%track is over
                %compute threshold;
                thresh = median(snr_revs(end-cur_revs-1:end)) + No
                break
            end
            cur_revs = 0;
            i_step_rev = i_step_rev + 1;
            
        end
        
    end
    thresh_mat(i_track) = thresh; %save the threshold for each of the tracks
end
median_thresh = median(thresh_mat) %compute the median threshold across all tracks in ntrack


%% Plot details from last track
figure
plot(t,template_ndiff)
xlabel('Time [sec]')
ylabel('model units')
title('Normalized Difference Template (from final track) see 1996b fig 12')
figure
hold on
plot(t,template_n)
plot(t,template_tin,'r')
xlabel('Time [sec]')
ylabel('model units')
legend('Noise-alone Template','Tone-plus-noise template')
title('Templates from final track')
figure
plot(t,noise_fadap_intn)
hold on
plot(t,noise3_fadap_intn,'k')
plot(t,tin_fadap_intn ,'r')
xlabel('Time [sec]')
ylabel('model units')
title('Internal representations from last trial')   
legend('Noise 1','Noise 2', 'Tone + Noise')
figure
subplot(2,1,1)
plot(t,template_ndiff)
legend('Normalized Difference Template')
title('Model inputs to decision device')
ylabel('model units')
subplot(2,1,2)
hold on
plot(t,noise_fadap_intn-template_n);
plot(t,noise3_fadap_intn-template_n,'k');
plot(t,tin_fadap_intn-template_n,'r');
legend('Noise 1','Noise 2', 'Tone + Noise')
xlabel('Time [sec]')
ylabel('model units')
title(['Each of the plots below is compared to the plot above Prob = ' num2str(prob)])
