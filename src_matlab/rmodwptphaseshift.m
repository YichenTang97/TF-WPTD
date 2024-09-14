function shiftedwpt = rmodwptphaseshift(wpt,Lo,Hi,level)
%   Adopted from and reversed the modwptphaseshift function defined under modwpt.m (line 372-438)

Lo = Lo/sqrt(2);
Hi = Hi/sqrt(2);
Lo = Lo(:)';
Hi = Hi(:)';

%Determine the size of the wavelet packets
numnodes = 2^level;
levels = level;


%Determine the center of energy
L = numel(Lo);
eScaling = sum((0:L-1).*Lo.^2);
eScaling = eScaling/norm(Lo,2)^2;
eWavelet = sum((0:L-1).*Hi.^2);
eWavelet = eWavelet/norm(Hi,2)^2;


bitvaluehigh = zeros(1,numnodes);
bitvaluelow = zeros(1,numnodes);
shiftedwpt = zeros(size(wpt));

% Compute phase shifts
m = 1;
for jj = 1:numel(levels)
    J = levels(jj);
    for nn = 0:2^J-1
        bitvaluehigh(m) = bitReversal(J,nn);
        bitvaluelow(m) = 2^J-1-bitvaluehigh(m);
        m = m+1;
    end
end


pJN = round(bitvaluelow*eScaling+bitvaluehigh*eWavelet);

for nn = 1:numnodes
    shiftedwpt(nn,:) = circshift(wpt(nn,:),[0 pJN(nn)]);
end
end


function bitvalue = bitReversal(J,N)

L = J;
filtsequence = zeros(1,J);
while J>=1
    
    remainder = mod(N,4);
    if (remainder == 0 || remainder == 3)
        filtsequence(J) = 0;
        
    else
        filtsequence(J) = 1;
    end
    J = J-1;
    N = floor(N/2);
    
end

bitvalue = sum(filtsequence.*2.^(0:L-1));
end

