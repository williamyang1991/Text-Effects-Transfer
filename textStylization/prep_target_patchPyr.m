function trgPatchPyr = prep_target_patchPyr(imgPyr, uvPixSub, LvlInd, iLvl, optS)

% PREP_TARGET_PATCHPYR
%
% Prepare multi-scale target patches according to patch centers
%
% Input:
%   - imgPyr
%   - uvPixSub: patch centers
%   - LvlInd: desired scales
%   - iLvl: current level
%   - optS
% Output:
%   - trgPatchPyr

trgPatchPyr = cell(size(LvlInd,2), 1);
ind = 1;
for i = LvlInd
    ratio = size(imgPyr{i},1)/ size(imgPyr{iLvl},1);
    trgPatch = prep_target_patch(imgPyr{i}, uvPixSub*ratio,  optS);
    trgPatchPyr{ind} = trgPatch;
    ind = ind + 1;
end

end