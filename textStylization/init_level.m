function [validPix, uvPix] = init_level(trgdistCur, srcdistCur, prad, lockRandOn, lockRandSeed)

% INIT_LEVEL
% Functionality: 
%   get uvPixels, validPixels in the level
%
% uvPixels: locations of patches in the target image (p in the paper)
% validPixels: locations of patches in the source image  (q in the paper)

% Get  uvPixSub, uvPixInd
[imgH, imgW] = size(trgdistCur);

uvMaps = ones(imgH, imgW);
uvMaps([1:prad, end-prad+1:end], :) = 0;
uvMaps(:, [1:prad, end-prad+1:end]) = 0;
[rUv, cUv] = find(uvMaps);
uvPix.sub = cat(2, cUv, rUv)';
uvPix.ind = sub2ind([imgH, imgW], rUv, cUv)';
uvPix.mask = uvMaps;
uvPix.numUvPix = size(uvPix.ind, 2);
uvPix.dist = trgdistCur(uvPix.ind);

% Get validPixels

[imgH, imgW] = size(srcdistCur);

validMap = ones(imgH, imgW);
validMap([1:prad,end-prad+1:end], :) = 0;
validMap(:, [1:prad,end-prad+1:end]) = 0;
[rV, cV] = find(validMap);
validPix.sub = cat(2, cV, rV)';
validPix.ind = sub2ind([imgH, imgW], rV, cV)';
validPix.mask = validMap;
validPix.numValidPix = size(validPix.ind, 2);
if lockRandOn
    rng(lockRandSeed);
end
[validPix.dist, validPix.distInd] = sort(srcdistCur(validPix.ind)+rand(1,validPix.numValidPix)*0.1);

end