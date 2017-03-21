function img = voting(img, srcimgCur, NNF, uvPix, wDistPatch, wDistImg, optS)

% VOTING: update the image using NNF
%
% Input: 
%   - wDistPatch
%   - srcimgCur
%   - NNF
%   - uvPix
%   - wDistPatch
%   - wDistImg, optS
% Output:
%   - img

[imgH, imgW, nCh] = size(img);

% Prepare source patch
numUvPix = size(uvPix.ind, 2);
uvValid.ind = true(1, numUvPix);    uvValid.pos = 1:numUvPix;
srcPatch = prep_source_patch(srcimgCur, NNF.uvTform.data, optS);

% Apply bias correction
if(optS.useBiasCorrection)
    temp = NNF.uvBias.data;
    temp(:,NNF.uvPix.dist<=1.2) = 0;
    biasPatch = reshape(temp, 1, nCh, numUvPix);
    srcPatch = bsxfun(@plus, srcPatch, biasPatch);
end

% Patch weight
wDistPatchC = reshape(wDistPatch(optS.pMidPix, :), 1, 1, numUvPix);
srcPatch = bsxfun(@times, srcPatch, wDistPatchC);

% Compute weighted average from source patches
srcPatch = reshape(srcPatch, optS.pNumPix*nCh, numUvPix);
imgAcc = zeros(imgH, imgW, 3, 'single');
for i = 1:numUvPix
    imgAcc(NNF.trgPatchInd(:,i)) = imgAcc(NNF.trgPatchInd(:,i)) + srcPatch(:,i);
end
img = imgAcc./wDistImg(:,:,ones(1,1,nCh));

end

