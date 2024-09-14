function [] = simulate_EEG()
    %% Synthesizing EEG trials with (success) and without (failure) error-related potentials (ErrPs).
    % Synthesized trials will be saved in the format of simulated_EEG_SNR_[-20,-15,...,10].mat
    % Each file contains two variables: 'data' which contains only ERPs and 'noisy_data' with noise added.
    % Both variables have shapes (300 trials x 32 channels x 307 timepoints). 
    % The first 150 trials correspond to success trials (no ErrPs) and the last 150 trials correspond to 
    % failure trials (with ErrPs).

    load('estimated_psd.mat');
    load('spatial_mapping.mat');

    % General parameters of the signal
    FRAMES = 307;
    TRIALS = 150; % num trials per class (succ/fail)
    SRATE = 256;

    %% Parameters for success trials
    % Parameters for N1
    N1AMP = 4;
    N1FREQ = 8;
    N1POS = 74; % 90 ms after stimuli onset

    % Parameters for P2
    P2AMP = 6;
    P2FREQ = 4;
    P2POS = 97; % 180 ms after stimuli onset

    % Parameters for P3
    P3AMP = 5;
    P3FREQ = 2;
    P3POS = 128; % 300 ms after stimuli onset

    %% Parameters for failure trials
    % Parameters for ERN
    ERNAMP = 4;
    ERNFREQ = 4;
    ERNPOS = 97; % 180 ms after stimuli onset

    % Parameters for P3a
    P3aAMP = 4;
    P3aFREQ = 2;
    P3aPOS = 128; % 300 ms after stimuli onset

    %%temporal jitter of peaks and amplitude of noise
    TJITTER_N1 = 3; % 12 ms
    TJITTER_P2 = 6; % 24 ms
    TJITTER_P3 = 9; % 36 ms

    TJITTER_ERN = 6; % 24 ms
    TJITTER_P3a = 9; % 36 ms

    for SNR=[-20,-15,-10,-5,0,5,10]
        fprintf('Generating data for SNR=%d dB...\n', SNR);

        disp ('Generating success trials');
        n1 = N1AMP * peak (FRAMES, TRIALS, SRATE, N1FREQ, N1POS, TJITTER_N1);
        p2 = P2AMP * peak (FRAMES, TRIALS, SRATE, P2FREQ, P2POS, TJITTER_P2);
        p3 = P3AMP * peak (FRAMES, TRIALS, SRATE, P3FREQ, P3POS, TJITTER_P3);
        data_succ = succ_mapping(:, 1) * n1 ...
            + succ_mapping(:, 2) * p2 ...
            + succ_mapping(:, 3) * p3;
        
        disp ('Generating failure trials')
        n1 = N1AMP * peak (FRAMES, TRIALS, SRATE, N1FREQ, N1POS, TJITTER_N1);
        p2 = P2AMP * peak (FRAMES, TRIALS, SRATE, P2FREQ, P2POS, TJITTER_P2);
        p3 = P3AMP * peak (FRAMES, TRIALS, SRATE, P3FREQ, P3POS, TJITTER_P3);
        ern = ERNAMP * peak (FRAMES, TRIALS, SRATE, ERNFREQ, ERNPOS, TJITTER_ERN);
        p3a = P3aAMP * peak (FRAMES, TRIALS, SRATE, P3aFREQ, P3aPOS, TJITTER_P3a);
        data_fail = succ_mapping(:, 1) * n1 ...
            + succ_mapping(:, 2) * p2 ...
            + succ_mapping(:, 3) * p3 ...
            + errp_mapping(:, 1) * ern ...
            + errp_mapping(:, 2) * p3a;
        
        data = [data_succ data_fail];
        
        disp ('Generating noise')
        noise_data = zeros([size(data)]);
        for ch = 1:size(data,1)
        noise_data(ch,:) = noise_data(ch,:) ...
            + gen_noise (FRAMES, TRIALS*2, SRATE, epsds(ch, :));
        end
        
        % Scale noise to match targeting SNR
        noise_data = (noise_data/std(noise_data,0,"all"))*std(data,0,"all")/db2mag(SNR);
        
        true_SNR = snr(data, noise_data);
        disp ('Actual SNR=');
        disp (true_SNR);
        
        noisy_data = data + noise_data;
        
        data = permute(reshape(data, [size(data,1),FRAMES,TRIALS*2]), [3,1,2]);
        noisy_data = permute(reshape(noisy_data, [size(noisy_data,1),FRAMES,TRIALS*2]), [3,1,2]);
        
        save(sprintf('simulated_EEG_SNR_%d.mat', round(SNR)), 'data', 'noisy_data');
        
        disp ('Done');
    end
end
