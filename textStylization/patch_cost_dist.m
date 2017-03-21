function [costApp, weight] = patch_cost_dist(trgPatch, srcPatch, optS)

% PATCH_COST_DIST
%
% Compute the weighted sum of the squared difference between
% cost between source and target patches

% Input:
%   - trgPatch, srcPatch, optS
% Output:
%   - costApp
%   - weight


% Initialization

numUvValidPix = size(trgPatch, 3);

% only the distance of the patch center is used
patchDist = trgPatch(optS.pMidPix,:,:) - srcPatch(optS.pMidPix,:,:);

% Sum of squared distance (dist(p)-dist(q))^2
patchDist = patchDist.^2;

% dist(p)*dist(p)
weight = reshape(trgPatch(optS.pMidPix,:,:), 1, numUvValidPix);

% distribution term cost (Eq. (13) in the paper)
costApp = reshape(patchDist, 1, numUvValidPix)./max(weight.^2,1);

end