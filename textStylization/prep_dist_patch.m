function [wDistPatch, wDistImg] = prep_dist_patch(distMap, trgPixPos, iLvl, optS)

% PREP_DIST_PATCH
% 
% Precompute patch weights and summation for patch matching and voting
% 
% Input:
%   - distMap
%   - trgPixPos
%   - iLvl
%   - optS
% Output:
%   - wDistPatch
%   - wDistImg

[imgH, imgW] = size(distMap);

wDistPatch  = prep_target_patch(distMap, trgPixPos, optS);
wDistPatch = bsxfun(@minus, wDistPatch, wDistPatch(optS.pMidPix,:));
wDistPatch  = optS.wDist(iLvl).^ (- wDistPatch); 

numUvPix = size(wDistPatch, 2);

wDistImg = zeros(imgH, imgW, 'single');
indMap = reshape(1:imgH*imgW, imgH, imgW);
indPatch  = prep_target_patch(indMap, trgPixPos, optS);

for i = 1: numUvPix
    wDistImg(indPatch(:,i)) = wDistImg(indPatch(:,i)) + wDistPatch(optS.pMidPix,i);
end
 
end