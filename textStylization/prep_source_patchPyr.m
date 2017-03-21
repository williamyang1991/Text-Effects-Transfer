function srcPatchPyr = prep_source_patchPyr(imgPyr, uvTform, LvlInd, iLvl, optS)

% PREP_SOURCE_PATCHPYR
%
% Prepare multi-scale source patches according to uvTform
%
% Input:
%   - imgPyr
%   - uvTform 
%   - LvlInd: desired scales
%   - iLvl: current level
%   - optS
% Output:
%   - srcPatchPyr


srcPatchPyr = cell(size(LvlInd,2), 1);
ind = 1;
uvTformCur = uvTform;
for i = LvlInd
    ratio = size(imgPyr{i},1)/ size(imgPyr{iLvl},1);
    uvTformCur(1,:) = uvTform(1,:)*ratio;
    uvTformCur(2,:) = uvTform(2,:)*ratio;
    srcPatch = prep_source_patch(imgPyr{i}, uvTformCur, optS);
    srcPatchPyr{ind} = srcPatch;
    ind = ind + 1;
end

end