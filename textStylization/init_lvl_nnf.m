function     [trgimgPyr, NNF, wDistPatch, wDistImg] = init_lvl_nnf(trgimgPyr, srcimgCur, NNF, mask, trgdistCur, srcdistCur, pSizehist, iLvl, optS)

% INIT_LVL_NNF
%
% Initialize the nearest neighbor field for the current level
%
% Input: 
%   - trgimgPyr, srcimgCur, NNF, mask, trgdistCur, srcdistCur, pSizehist, iLvl, optS
% Output:
%   - trgimgPyr
%   - NNF 
%   - wDistPatch
%   - wDistImg

% Prepare distance weight
mask(mask < 0.3) = 0;
[distMap, ~] = bwdist(mask, 'euclidean');

if(iLvl == optS.numPyrLvl)
    % Initialize the NNF for the coarest level using distance-based random sampling 
    NNF = init_nnf(trgdistCur, srcdistCur, pSizehist, optS);
    % patch weighting
    [wDistPatch, wDistImg] = prep_dist_patch(distMap, NNF.uvPix.sub, iLvl, optS);
    % update the image using NNF
    trgimgPyr{iLvl} = voting(trgimgPyr{iLvl}, srcimgCur, NNF, NNF.uvPix, wDistPatch, wDistImg, optS);
else
    % Initialize the NNF upsampling of NNF from previous level
    NNF = upsample(NNF, trgdistCur, srcdistCur, pSizehist, optS);
    [wDistPatch, wDistImg] = prep_dist_patch(distMap, NNF.uvPix.sub, iLvl, optS);
    trgimgPyr{iLvl} = voting(trgimgPyr{iLvl}, srcimgCur, NNF, NNF.uvPix, wDistPatch, wDistImg, optS);
end

% get multi-scale images (used for the Appearance Term)
for i = optS.numPyrLvl:-1:iLvl+1
    trgimgPyr{i} = imresize(trgimgPyr{iLvl}, [size(trgimgPyr{i},1), size(trgimgPyr{i},2)], optS.resampleKernel);
end

NNF.uvDtBdPixPos = double(distMap(NNF.uvPix.ind));

end