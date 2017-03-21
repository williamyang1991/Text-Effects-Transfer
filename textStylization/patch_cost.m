function [costPatchCand, uvBiasCand] = ...
    patch_cost(trgPatchPyr, srcPatchPyr, trgTextPatchPyr, srcTextPatchPyr, trgDistPatch, srcDistPatch, wDistPatchCur, ...
    freqMap, srcInd, pSizeweight, optS, iLvl, iter, numIterLvl)

% PATCH_COST
% 
% compute the objective function (Eq. (10) in the paper)
%
% Input:
%   - trgPatchPyr, trgTextPatchPyr, trgDistPatch: P(p), P'(p), dist(p)
%   - srcPatchPyr, srcTextPatchPyr, srcDistPatch: Q(q), Q'(q), dist(q) 
%   - wDistPatchCur: voting weight for patches
%   - NNF    
%   - wDistPatch
%   - optS
%   - iLvl, numIterLvl: current level, current iteration
%   - lockAngleFlag: whether to use patch rotation
% Output:
%   - NNF
%   - trgimgPyr


% srcInd 指示每个目标块的源块的索引

numUvPix = size(wDistPatchCur, 2);
costPatchCand = zeros(4, numUvPix);
[costApp, uvBiasCand] = patch_cost_app(trgPatchPyr{1}, srcPatchPyr{1}, wDistPatchCur, optS, optS.useBiasCorrection);
costApp = costApp .* pSizeweight(1,:);
costAppText = patch_cost_app(trgTextPatchPyr{1}, srcTextPatchPyr{1}, wDistPatchCur, optS);
costAppText = costAppText .* pSizeweight(1,:);

% Patch cost - appearance
for i = 2:size(trgPatchPyr,1)    
    [tmpcostApp, ~] = patch_cost_app(trgPatchPyr{i}, srcPatchPyr{i}, wDistPatchCur, optS, optS.useBiasCorrection, uvBiasCand);    
    costApp = costApp + tmpcostApp.*pSizeweight(i,:);
    [tmpcostAppText, ~] = patch_cost_app(trgTextPatchPyr{i}, srcTextPatchPyr{i}, wDistPatchCur, optS);
    costAppText = costAppText + tmpcostAppText.*pSizeweight(i,:);
end

% Patch cost - distribution
[costAppDist, distweight] = patch_cost_dist(trgDistPatch, srcDistPatch, optS);

% Patch cost - repetitiveness
freqData = freqMap(:)';
costRep = freqData(srcInd);

% Weighted sum of the costs
costPatchCand(1,:) = costApp;
costPatchCand(2,:) = optS.lambdaText * costAppText;
costPatchCand(3,:) = optS.lambdaDist * costAppDist;
costPatchCand(4,:) = optS.lambdaRep * costRep; 

% Adaptive weighting for patch cost
% loose the penalty of character shape differences along with the iteration
costPatchCand(2,:) = (1-iter/numIterLvl)*costPatchCand(2,:);
% clamp distribution term
costPatchCand(3,:) = min(1, costPatchCand(3,:));
% normalize repetitiveness based on the size of S and T
% loose the penalty of repetitiveness as the distance to the text goes further
costPatchCand(4,:) = optS.repCostRatio * costPatchCand(4,:) ./ max(1,abs(distweight));

end
