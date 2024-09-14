function Xdwpt = modwptdec(X, LoD, HiD, level)
%MODWPTDEC Decompose signals.
%   X with shape n_trials x n_channels x n_times

    nTrials = size(X,1);
    nChs = size(X,2);
    nTimes = size(X,3);
    nFreqs = 2^level;

    Xdwpt = zeros(nTrials, nChs, nFreqs, nTimes);

    for x = 1:nTrials
        for ch = 1:nChs
            Xdwpt(x,ch,:,:) = modwpt(squeeze(X(x,ch,:))', LoD, HiD, level, 'TimeAlign', true);
        end
    end
end

