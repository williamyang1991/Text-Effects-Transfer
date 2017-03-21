function [trgimgPyr, imgPyrNNF] = ...
    synthesis(srcimgPyr, srctextPyr, srcdistPyr, trgimgPyr, trgtextPyr, trgdistPyr, pSizehist, optS)

% SYNTHESIS:
%
% Patch-based synthesis using patchmatch algorithm
%
% Input:
%   - srcimgPyr, srctextPyr, srcdistPyr, trgimgPyr, trgtextPyr, trgdistPyr
%   - pSizehist
%   - optS
% Output:
%   - trgimgPyr
%   - imgPyrNNF
%

imgPyrNNF = cell(optS.numPyrLvl, 1);
% nearest neighbor field, for each p, store the matched q
NNF = [];
numIterLvl = optS.numIter;
pyrLvl = optS.numPyrLvl: -1 : optS.topLevel;

% Coarse-to-fine text effects transfer
for iLvl = pyrLvl
    % Initialize level   
    srcdistCur = srcdistPyr{iLvl};
    trgdistCur   = trgdistPyr{iLvl};    
    srcimgCur = srcimgPyr{iLvl};
    mask = trgtextPyr{iLvl};
    
    % === Prepare img and NNF for the current level ===
    fprintf('--- Initialize NNF: ');
    [trgimgPyr, NNF, wDistPatch, ~] = ...
        init_lvl_nnf(trgimgPyr, srcimgCur, NNF, mask, trgdistCur, srcdistCur, pSizehist, iLvl, optS);
    
    % normalize the Psycho-Visual Term 
    % used when there is great size difference for S and T
    optS.repCostRatio = (NNF.validPix.numValidPix / NNF.uvPix.numUvPix).^2;
      
    % Number of iterations at the currect level
    numIterLvl = max(numIterLvl - optS.numIterDec, optS.numIterMin);
    
    fprintf('--- Pass... level: %d, #Iter: %d, #uvPixels: %7d\n', iLvl, numIterLvl, NNF.uvPix.numUvPix);
    fprintf('--- %3s\t%12s\t%12s\t%12s\t%10s\n', 'iter', '#PropUpdate', '#RandUpdate', '#RegUpdate', 'AvgCost');
    
    % patchmatch (interation of random search and propagation)
    if(iLvl == optS.numPyrLvl)
        % for the first iteration, patch roration is not used
        [trgimgPyr, NNF] = one_pass(trgimgPyr, trgtextPyr, trgdistCur, srcimgPyr, srctextPyr, srcdistCur, NNF, wDistPatch, numIterLvl, iLvl, optS, 1);
        % for the following iterations, patch roration is used
        [trgimgPyr, NNF] = one_pass(trgimgPyr, trgtextPyr, trgdistCur, srcimgPyr, srctextPyr, srcdistCur, NNF, wDistPatch, numIterLvl, iLvl, optS, 0);         
    else
        [trgimgPyr, NNF] = one_pass(trgimgPyr, trgtextPyr, trgdistCur, srcimgPyr, srctextPyr, srcdistCur, NNF, wDistPatch, numIterLvl, iLvl, optS, 0);
    end
    
    fprintf('Max cost & min cost: %f, %f\n', max(NNF.uvCost.data), min(NNF.uvCost.data));
    % Save the result
    imgPyrNNF{iLvl} = NNF;
end

end