function [den_coeffs] = NZT_fit(x, stim, scale)
%%% fitting NZT
% params:
% x: EEG data in trials x channels x times
% stim: the sample index for the stimuli onset
% scale: levels of wavelet decomposition to perform

n_chans = size(x,2);
n_times = size(x,3);

% zero padding
expected_len = 2^scale * ceil(n_times/(2^scale));
if n_times < expected_len
    x(:,:,n_times+1:expected_len) = 0;
end

% denoise by channel
den_coeffs = cell(1,n_chans);
for i=1:n_chans
    av = squeeze(mean(x(:,i,:),1));
    [~,~,den_coeff,~,~]=Run_NZT(av,stim,scale);
    den_coeffs{i} = den_coeff;
end
end