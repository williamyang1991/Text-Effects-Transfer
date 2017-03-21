function uvPix = get_pixel_ind(trgdistCur, prad)

% GET_PIXEL_IND
% Functionality: 
%   prepare patch centers

[imgH, imgW, ~] = size(trgdistCur);

% boundary treatment
uvMaps = ones(imgH, imgW);
uvMaps([1:prad, end-prad+1:end], :) = 0;
uvMaps(:, [1:prad, end-prad+1:end]) = 0;
[rUv, cUv] = find(uvMaps);
% patch centers
uvPix.sub = cat(2, cUv, rUv)';
uvPix.ind = sub2ind([imgH, imgW], rUv, cUv)';
uvPix.mask = uvMaps;
uvPix.numUvPix = size(uvPix.ind, 2);

end