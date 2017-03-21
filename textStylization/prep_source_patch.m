function srcPatch = prep_source_patch(img, uvTform, optS)

% PREP_SOURCE_PATCH
%
% Prepare source patches according to uvTform
%
% Input:
%   - img
%   - uvTform 
%   - optS
% Output:
%   - srcPatch


numUvPix = size(uvTform, 2);

srcPatchPos = zeros(optS.pNumPix, 2, numUvPix);

% rotation
sintheta = sin(pi * uvTform(3,:));
costheta = cos(pi * uvTform(3,:));

% pixel potisions of the source patches
for i = 1 : optS.pNumPix
    dx = optS.refPatchPos(1,i);
    dy = optS.refPatchPos(2,i);    
    srcPatchPos(i, 1, :) = costheta*dx - sintheta*dy + uvTform(1,:);
    srcPatchPos(i, 2, :) = sintheta*dx + costheta*dy + uvTform(2,:);
end

% Avoid sample out of boundary positions
srcPatchPos(:,1,:) = clamp(srcPatchPos(:,1,:), 1, size(img,2));
srcPatchPos(:,2,:) = clamp(srcPatchPos(:,2,:), 1, size(img,1));

% Sample target patch
img = im2double(img);
srcPatch = mirt2D_mexinterp(img, srcPatchPos(:, 1, :), srcPatchPos(:, 2, :));
if(size(srcPatch, 3) > 1)
    srcPatch = permute(srcPatch, [1, 3, 2]);
end

end