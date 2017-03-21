function [imgPyr, scaleImgPyr] = create_img_pyramid(img, optS)

% CREATE_IMG_PYRAMID
%
% Create image pyramid with linear or log scale for coarse to fine image
% completion
%
% Input:
%   - img:  Image
%   - optS: options
% Output:
%   - imgPyr:      Image pyramid
%   - scaleImgPyr: Image dimensions in each level

% Image size in the high-resolution image
[imgHeight, imgWidth, nCh] = size(img);

% Compute the coarsest image scale
imgSizeMin = min(imgHeight, imgWidth);
coarestScale = optS.coarestImgSize/imgSizeMin;

% Compute the scale in each layer in the image pyramid
if(optS.useLogScale) % use log scale
    scalePyr = 2.^linspace(0, log2(coarestScale), optS.numPyrLvl);
else % use linear scale
    scalePyr = linspace(1, coarestScale, optS.numPyrLvl);
end

% Image size in each layer
imgHPyr = round(imgHeight *scalePyr);
imgWPyr = round(imgWidth  *scalePyr);

% Initialize image pyramid
imgPyr  = cell(optS.numPyrLvl, 1);
scaleImgPyr = cell(optS.numPyrLvl, 1);

% Finest level (high-resolution image)
imgPyr{1} = img;
scaleImgPyr{1}.imgScale = 1;
scaleImgPyr{1}.imgSize = [imgHeight, imgWidth];

% Downsampled images
for k = 2: optS.numPyrLvl
    imgPyr{k} = imresize(img, [imgHPyr(k), imgWPyr(k)], optS.resampleKernel);
    scaleImgPyr{k}.imgScale = scalePyr(k);
    scaleImgPyr{k}.imgSize  = [imgHPyr(k), imgWPyr(k)];
end

end