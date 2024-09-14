function denoisedX = modwptdenoise(X, selectors, LoD, HiD, level)
%MODWPTDENOISE Decompose, select significant coefficients, and reconstruct.
%   X in the shape of n_trials x n_chs x n_times
%   selectors should be a cell array of shape 1 x n_selectors. Selectors 
%       must be the index for the flattened Xdwpt (decomposed X) in 
%       Fortran-like index order. Each selector is a vector.

    nTrials = size(X,1);
    nChs = size(X,2);
    nTimes = size(X,3);
    nFreqs = 2^level;

    Xdwpt = modwptdec(X,LoD,HiD,level);

    % assume multiple selectors been used if selectors is a cell array

    denoisedX = zeros(nTrials, nChs*size(selectors,2), nTimes);
    for i = 1:size(selectors,2)
        tempX = zeros(size(Xdwpt));
        for x = 1:nTrials
            xdwptFlat = reshape(Xdwpt(x,:,:,:), 1, []);
            xdwptSel = zeros(size(xdwptFlat));
            xdwptSel(selectors{i}) = xdwptFlat(selectors{i});
            tempX(x,:,:,:) = reshape(xdwptSel,nChs,nFreqs,nTimes);
        end
        denoisedX(:,(i-1)*nChs+1:i*nChs,:) = modwptrec(tempX,LoD,HiD,level);
    end
end

