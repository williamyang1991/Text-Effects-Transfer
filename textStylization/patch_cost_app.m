function [costApp, uvBias] = patch_cost_app(trgPatch, srcPatch, wDistPatch, optS, useBiasCorrection, uvBias)

% PATCH_COST_APP
%
% Compute the weighted sum of the squared difference between
% cost between source and target patches

% Input:
%   - trgPatch, srcPatch, wDistPatch, optS, useBiasCorrection, uvBias
% Output:
%   - costApp
%   - uvBias

% Initialization

if nargin == 4
    useBiasCorrection = false;
    uvBias = zeros(3, size(wDistPatch, 2));
elseif nargin == 5
    uvBias = zeros(3, size(wDistPatch, 2));
end

numUvValidPix = size(wDistPatch, 2);

% Apply bias correction
if(useBiasCorrection)
    % Mean of source and target patch
    if nargin == 5
        meanTrgPatch = mean(trgPatch, 1);
        meanSrcPatch = mean(srcPatch, 1);

        % Compute bias and clamp it to inteval [optS.minBias, optS.maxBias]
        biasPatch = meanTrgPatch - meanSrcPatch;
        biasPatch = clamp(biasPatch, optS.minBias, optS.maxBias);

        % Update the UV map for gain and bias
        uvBias = reshape(biasPatch, 3, numUvValidPix);
    else
        biasPatch = reshape(uvBias, 1, 3, numUvValidPix);
    end
    % Compute patch appearance cost
    srcPatch = bsxfun(@plus, srcPatch, biasPatch);
end
patchDist = trgPatch - srcPatch;

% Weight patch appearance cost by their distances to text
wDistPatchC = reshape(wDistPatch, optS.pNumPix, 1, numUvValidPix);

% Sum of squared distance
if(strcmp(optS.costType, 'L1'))
    patchDist = abs(patchDist);
elseif(strcmp(optS.costType, 'L2'))
    patchDist = patchDist.^2;
end

% Apply weights
patchDist = bsxfun(@times, patchDist, wDistPatchC);
patchDist = sum(sum(patchDist, 1),2);
patchDist = reshape(patchDist, 1, numUvValidPix);

% Weight normalization
sumDistWeight = sum(wDistPatch, 1);%1*N
costApp = patchDist./sumDistWeight;%1*N

end