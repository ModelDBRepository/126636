function nout = noise_gen(fs,dur,lofreq,hifreq,No,avrms,rampdur);
%Creates a bandpass noise with infinite rejection in the frequency domian
% 
% nout = noise_gen(fs,dur,lofreq,hifreq,No,avrms,rampdur);
% 
% No, the spectrum level;avrms, the average rms;
% and rampdur, the ramp duration in seconds, are optional.
    tdres = 1/fs;
    totalpts = round(dur.*fs);
  
    stim1 = randn(1,totalpts);  
    
    fres = 1./(totalpts * tdres);
    
    % Now filter the noise 
    stimspectrum = fft(stim1);
    nptslofreq = floor(lofreq/fres);
    nptshifreq = ceil(hifreq/fres);
    stimspectrum(1:nptslofreq) = 0.;  % zero out all bands that apply 
    stimspectrum(nptshifreq:(totalpts - nptshifreq+2)) = 0.;
    stimspectrum((totalpts - nptslofreq+2):totalpts) = 0.;
    % So, values in the fft of stim1 (both real & imag parts) have been zeroed
    out = real(ifft(stimspectrum)); % inverse FFT to create noise waveform
    
    if nargin == 5;
        out = out/rms(out);
        nout = out * 20e-6* 10^((No + 10*log10(hifreq-lofreq))/20);
    elseif nargin > 5
        out = out/avrms;
        nout = out * 20e-6* 10^((No + 10*log10(hifreq-lofreq))/20);
      if nargin ==7
        rmp = hanning(round(rampdur*2*fs));
        nout = nout .* [rmp(1:round(length(rmp)/2))' ones(1,length(nout)-length(rmp)) rmp(round(length(rmp)/2)+1:end)'];    
      end
    else 
        nout = out;
    end