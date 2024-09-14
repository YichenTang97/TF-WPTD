function [x_filt] = NZT_transform(x,den_coeffs,scale)
%%% filtering using fitted coefficients
% params:
% x: EEG data in trials x channels x times
% den_coeffs: the denoising coefficients found using NZT_fit
% scale: levels of wavelet decomposition to perform, should be the same
% with the one used in NZT_fit

[n_samps,n_chans,n_times] = size(x);

% zero padding
expected_len = 2^scale * ceil(n_times/2^scale);
if n_times < expected_len
    x(:,:,n_times+1:expected_len) = 0;
end

x_filt = zeros([n_samps,n_chans,n_times]);
for i=1:n_chans
    x_chan = squeeze(x(:,i,:));
    den_coeff = den_coeffs{i};
    x_chan_filt = st_den(reshape(x_chan.', [numel(x_chan),1]), den_coeff, scale, expected_len);
    x_filt(:,i,1:n_times) = x_chan_filt(:,1:n_times);
end
end

