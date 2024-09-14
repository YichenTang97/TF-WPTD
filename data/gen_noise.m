function signal = gen_noise (frames, epochs, srate, epsd)

% function signal = gen_noise (frames, epochs, srate, epsd)
%
% Function generates noise with the power spectrum of human EEG
% Inputs:
%  frames - number of signal frames per each trial
%  epochs - number of simulated trials
%  srate - sampling rate of simulated signal
%  epsd - power spectrum density of EEG noise
% Output:
%  signal - simulated EEG signal; vector: 1 by frames*epochs containing concatenated trials
% Edited: XXXX XXXX, Oct, 2022
% Adopted from the original noise function implemented by: Rafal Bogacz and Nick Yeung, Princeton Univesity, December 2002

sumsig = 200;	%number of sinusoids from which each simulated signal is composed of

signal = zeros (1, epochs * frames);
for trial = 1:epochs
   freq=0;
   range = [(trial-1)*frames+1:trial*frames];
   for i = 1:sumsig
      freq = freq + (1*rand(1));
      freqamp = epsd(min (ceil(freq), 128)) / epsd(1);
      phase = rand(1)*2*pi;
      signal (range) = signal (range) + sin ([1:frames]/srate*2*pi*freq + phase) * freqamp;
   end
end