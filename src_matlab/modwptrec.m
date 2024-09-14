function X = modwptrec(Xdwpt, LoD, HiD, level)
%MODWPTREC Reconstruct Xdwpt.
%   Xdwpt in the shape of n_trials x n_channels x n_freqs x n_times

    nTrials = size(Xdwpt,1);
    nChs = size(Xdwpt,2);
    nFreqs = size(Xdwpt,3);
    nTimes = size(Xdwpt,4);
    
    X = zeros(nTrials, nChs, nTimes);

    for x = 1:nTrials
        for ch = 1:nChs
            unshiftedXDwpt = rmodwptphaseshift(squeeze(Xdwpt(x,ch,:,:)),LoD,HiD,level);
            X(x,ch,:) = imodwpt(unshiftedXDwpt, LoD, HiD);
        end
    end
end

