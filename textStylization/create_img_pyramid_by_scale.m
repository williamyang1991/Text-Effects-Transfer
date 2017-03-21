function imgPyr = create_img_pyramid_by_scale(img, optS, scaleImgPyr)

% CREATE_IMG_PYRAMID_BY_SCALE
%
% Create image pyramid according to the given downsampling scales
%
% Input:
%   - img:  Image
%   - optS: options
%   - scaleImgPyr: given downsampling scales
% Output:
%   - imgPyr:      Image pyramid

% Initialize image pyramid
imgPyr  = cell(optS.numPyrLvl, 1);

% Finest level (high-resolution image)
imgPyr{1} = img;

% Downsampled images
for k = 2: optS.numPyrLvl   
    imgPyr{k} = imresize(img, scaleImgPyr{k}.imgScale, optS.resampleKernel);
end

end