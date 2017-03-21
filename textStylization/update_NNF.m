function [NNF, nUpdate]= update_NNF(trgPatchPyr, srcimgPyr, trgTextPatchPyr, srcTextPyr, trgDistPatch, srcDist, wDistPatch, NNF, iLvl, optS, iter, numIterLvl, lockAngleFlag)

% UPDATE_NNF
%
% Update the nearest neighbor field using the PatchMatch algorithm
%
% Input:
%   - trgPatchPyr, trgTextPatchPyr, trgDistPatch: P(p), P'(p), dist(p)
%   - srcimgPyr, srcTextPyr, srcDist: source of Q(q), Q'(q), dist(q) 
%   - wDistPatch: voting weight for patches
%   - NNF       
%   - optS
%   - iLvl, iter, numIterLvl: adjust cost weights based on level and iteration
%   - lockAngleFlag: whether to use patch rotation
% Output:
%   - NNF
%   - nUpdate

nUpdate = zeros(1,3);

for i = 1:optS.numPassPerIter
    % propagate along four directions
    for iDirect = 1:4
        [NNF, n] = propagate(trgPatchPyr, srcimgPyr, trgTextPatchPyr, srcTextPyr, trgDistPatch, srcDist, ...
                    wDistPatch, NNF, optS, iDirect, iLvl, iter, numIterLvl);
        nUpdate(1) = nUpdate(1) + n;
    end
    
    if(iLvl > optS.propOnlyLevel)        
        % Random sampling
        [NNF, n] = random_search(trgPatchPyr, srcimgPyr, trgTextPatchPyr, srcTextPyr, trgDistPatch, srcDist, ...
                    wDistPatch, NNF, optS, iLvl, iter, numIterLvl, lockAngleFlag);
        nUpdate(2) = nUpdate(2) + n;
    end
end

end