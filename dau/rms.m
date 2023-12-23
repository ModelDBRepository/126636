function [rms_val] = rms(wave)
rms_val = sqrt(mean(wave .^ 2));