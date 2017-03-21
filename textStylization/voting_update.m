function img = voting_update(wDistPatch, img, srcimgCur, NNF, optS, useBiasCorrection)

% VOTING_UPDATE
% 
% Fast voting update, instead of computing the weighted average from all
% source patches, here we only update patches that have been updated
%

if nargin == 5
    useBiasCorrection = false;
end

[imgH, imgW, nCh] = size(img);
if nCh == 1
    img = img(:,:,ones(3,1));
    srcimgCur = srcimgCur(:,:,ones(3,1));
    nCh = 3;
end
numUvPix = sum(NNF.update.data);
if(numUvPix~=0)
    
    % Prepare source patch
    uvTform = NNF.uvTform.data(:, NNF.update.data);
    srcPatch = prep_source_patch(srcimgCur, uvTform, optS);
    
    if(useBiasCorrection)
        biasPatch = NNF.uvBias.data(:, NNF.update.data);
    end
    
    wDistPatch = wDistPatch(:, NNF.update.data);
    trgPatchInd = NNF.trgPatchInd(:, NNF.update.data);
    
    % Apply bias correction
    if(useBiasCorrection)
        biasPatch = reshape(biasPatch, 1, nCh, numUvPix);
        srcPatch = bsxfun(@plus, srcPatch, biasPatch);
    end
    
    % Patch weight
    wDistPatchC = reshape(wDistPatch(optS.pMidPix, :), 1, 1, numUvPix);
    srcPatch = bsxfun(@times, srcPatch, wDistPatchC);
        
    % Compute weighted average from source patches
    srcPatch = reshape(srcPatch, optS.pNumPix*nCh, numUvPix);
    imgAcc = zeros(imgH, imgW, 3, 'single');
    weightAcc = optS.voteUpdateW*ones(imgH, imgW, 'single');
    
    for i = 1:numUvPix
        imgAcc(trgPatchInd(:,i)) = imgAcc(trgPatchInd(:,i)) + srcPatch(:,i);
        weightAcc(trgPatchInd(1:optS.pNumPix,i)) = weightAcc(trgPatchInd(1:optS.pNumPix,i)) + wDistPatch(optS.pMidPix, i);
    end
    imgAcc = imgAcc + optS.voteUpdateW*img;
    img = imgAcc./weightAcc(:,:,ones(1,1,nCh));
end

end