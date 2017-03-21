function NNFVis= vis_nnf(NNF)

% VIS_NNF
%     visualization of NNF
% 

mask = NNF.uvPix.mask;

[imgH, imgW] = size(mask);

bdMask = edge(mask);
bdMask = bdMask(:,:,ones(3,1));
mask = mask(:,:,ones(3,1)); 

%% uvTfomMapVis
NNFVis.uvTfomMapVis = vis_tform_map(NNF.uvTform.map, mask, bdMask);

%% angleMapVis

NNFVis.angleMapVis = 0.5+max(min(0.5, NNF.uvTform.map(:,:,3)),-0.5);
NNFVis.angleMapVis(1,1) = 0;
NNFVis.angleMapVis(end,end) = 1;

%% uvCostMapVis

NNFVis.uvCostMapVis = (NNF.uvCost.map - min(NNF.uvCost.map(:)))/max(NNF.uvCost.map(:));


%% uvPixUpdateSrc

NNFVis.uvPixUpdateSrcMap = zeros(imgH, imgW, 3);
for ch = 1:3
    NNFVis.uvPixUpdateSrcMap(:,:,ch) = im2double(NNF.uvPixUpdateSrc.map == ch);
end

end

function NNFMapVis = vis_tform_map(NNFMap, mask, bdMask)

[imgH, imgW, ch] = size(mask);

% Initialize NNFMapVis
NNFMapVis = zeros(imgH, imgW, 3, 'single');
[X, Y] = meshgrid(1:imgW, 1:imgH);
X = X/imgW;     Y = Y/imgH;

NNFMapVis(:,:,2) = 0.5;
NNFMapVis(:,:,1) = X;
NNFMapVis(:,:,3) = Y;

% Prepare the visualization of the NNF
NNFMapCurVis = zeros(imgH, imgW, 3, 'single');
NNFMapCurVis(:,:,2) = 0.5;
NNFMapCurVis(:,:,1) = NNFMap(:,:,1)/imgW;
NNFMapCurVis(:,:,3) = NNFMap(:,:,2)/imgH;

NNFMapVis = NNFMapVis.*(1-mask) + NNFMapCurVis.*mask;
NNFMapVis(bdMask) = 1;

end