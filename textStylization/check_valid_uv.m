function uvValidInd = check_valid_uv(srcPos, validSrcMask)

% CHECK_VALID_UV 
% Functionality: 
%   check patches whether lie in the valid region of validSrcMask
%


uvSub = round(srcPos);

% clamp
uvSub(1,:) = clamp(uvSub(1,:), 1, size(validSrcMask,2));
uvSub(2,:) = clamp(uvSub(2,:), 1, size(validSrcMask,1));

uvInd = sub2ind(size(validSrcMask), uvSub(2,:), uvSub(1,:));

uvValidInd = validSrcMask(uvInd);

end