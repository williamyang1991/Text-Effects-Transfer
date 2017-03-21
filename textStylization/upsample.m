function NNF_H = upsample(NNF_L, trgdistCur, srcdistCur, pSizehist, optS)
%
% UPSAMPLE: upsample the nearest neighbor field
%
% Input:
%     - NNF_L: nearest neighbor field of the low resolution image
%     - trgdistCur: distance map of T
%     - srcdistCur: distance map of S
%     - pSizehist: posterior probability
%     - optS: parameters for synthesis process
% Output:
%   Image height and width
%     - NNF.imgH: image height of T
%     - NNF.imgW: image width of T
%     - NNF.imgH2: image height of S
%     - NNF.imgW2: image width of S
%   Precomputing pixel positions for PatchMatch algorithm
%     - NNF.uvPix: the target patch positions p of T/T'
%     - NNF.uvPixN: the neighboring patch positions of 4 connected pixels
%     - NNF.validPix: the source patch positions q of S/S'
%     - NNF.trgPatchInd: the indice to retrieve target patches
%   Initialize NNF components:
%     - NNF.uvTform: the matched source patch positions for target patches
%         - NNF.uvTform.data:   3 x numUvPix
%         - NNF.uvTform.map:    H x W x 3
%     - NNF.uvBias: the bias between source and target patches
%         - NNF.uvBias.data:    3 x numUvPix
%         - NNF.uvBias.map:     H x W x 3
%     - NNF.uvCost: the matching cost between source and target patches
%         - NNF.uvCost.data:    1 x numUvPix
%         - NNF.uvCost.map:     H x W x 1
%     - NNF.freq: the number of pixels that find q as its correspondence
%         - NNF.freq.data:      1 x numValidPix
%         - NNF.freq.map:       H2 x W2 x 1
%     - NNF.freqCost: the cost of Psycho-Visual Term
%         - NNF.freqCost.data:  1 x numUvPix
%         - NNF.freqCost.map:   H x W x 1

%% === Initialize uvPix and validPix ===
[NNF_H.imgH, NNF_H.imgW] = size(trgdistCur);
[NNF_H.imgH2, NNF_H.imgW2] = size(srcdistCur);
[NNF_H.validPix, NNF_H.uvPix] = init_level(trgdistCur, srcdistCur, optS.pRad, optS.lockRandOn, optS.lockRandSeed);

%% === Initialize uvPixN ===
NNF_H.uvPixN = cell(4,1);
for i = 1: 4
    NNF_H.uvPixN{i}.sub(1,:) = NNF_H.uvPix.sub(1,:) - optS.propDir(1,i);
    NNF_H.uvPixN{i}.sub(2,:) = NNF_H.uvPix.sub(2,:) - optS.propDir(2,i);
    NNF_H.uvPixN{i}.ind = sub2ind([NNF_H.imgH, NNF_H.imgW], NNF_H.uvPixN{i}.sub(2,:), NNF_H.uvPixN{i}.sub(1,:));
    NNF_H.uvPixN{i}.validInd = NNF_H.uvPix.mask(NNF_H.uvPixN{i}.ind);
end

%% === Initialize indMap ===
NNF_H.pixIndMap = reshape(1:NNF_H.imgH*NNF_H.imgW, NNF_H.imgH, NNF_H.imgW);
indTrgPatch = prep_target_patch(NNF_H.pixIndMap, NNF_H.uvPix.sub, optS);

NNF_H.trgPatchInd = cat(1, indTrgPatch, ...
    indTrgPatch+NNF_H.imgH*NNF_H.imgW, indTrgPatch + 2*NNF_H.imgH*NNF_H.imgW);

%% === Initialize uvPixL ===

imgH_H = NNF_H.imgH;    imgW_H = NNF_H.imgW;
imgH_L = NNF_L.imgH;    imgW_L = NNF_L.imgW;

sX = imgH_L/imgH_H;     sY = imgW_L/imgW_H;
uvPixL.sub = round(diag([sX, sY])*NNF_H.uvPix.sub);
uvPixL.sub(1,:) = clamp(uvPixL.sub(1,:), optS.pRad+1, imgW_L - optS.pRad);
uvPixL.sub(2,:) = clamp(uvPixL.sub(2,:), optS.pRad+1, imgH_L - optS.pRad);
uvPixL.ind = sub2ind([imgH_L, imgW_L], uvPixL.sub(2,:), uvPixL.sub(1,:));


%% === Initialize uvPixUpdateSrc ===
NNF_H.uvPixUpdateSrc.map = zeros(NNF_H.imgH, NNF_H.imgW);
NNF_H.uvPixUpdateSrc.data = uvMat_from_uvMap(NNF_L.uvPixUpdateSrc.map, uvPixL);
NNF_H.uvPixUpdateSrc.map = update_uvMap(NNF_H.uvPixUpdateSrc.map, NNF_H.uvPixUpdateSrc.data, NNF_H.uvPix, true(1,NNF_H.uvPix.numUvPix));

%% === Update uvTform ===

uvTform_L = uvMat_from_uvMap(NNF_L.uvTform.map, uvPixL);
uvTform_L(1:2,:) = diag([1/sX, 1/sY])*uvTform_L(1:2,:);

% Refinement
refineVec = NNF_H.uvPix.sub - diag([1/sX, 1/sY])*uvPixL.sub;
uvTform_H = trans_tform(uvTform_L, refineVec);

% Clamp
uvTform_H(1,:) = clamp(uvTform_H(1,:), optS.pRad+1, imgW_H - optS.pRad);
uvTform_H(2,:) = clamp(uvTform_H(2,:), optS.pRad+1, imgH_H - optS.pRad);
uvValid_H = check_valid_uv(uvTform_H(1:2,:), NNF_H.validPix.mask);

uvInvalidInd = ~uvValid_H;
nInvalidUv_H = sum(uvInvalidInd);
if(nInvalidUv_H)
    randInd = randi(size(NNF_H.validPix.ind, 2), nInvalidUv_H, 1);
    uvRand = NNF_H.validPix.sub(:, randInd);
    uvTform_H(1:2,uvInvalidInd) = uvRand;
end

% Update uvTform.map
NNF_H.uvTform.data = uvTform_H;
I = zeros(1, 1, 3);
NNF_H.uvTform.map = repmat(I, [imgH_H, imgW_H, 1]);
NNF_H.uvTform.map = update_uvMap(NNF_H.uvTform.map, NNF_H.uvTform.data, NNF_H.uvPix, true(1,NNF_H.uvPix.numUvPix));

%% === Initialize source frequency ===

[NNF_H.freq, NNF_H.freqCost] = get_NNF_freq(NNF_H);
NNF_H.freqCost.data = NNF_H.freqCost.data;
NNF_H.freqCost.map = NNF_H.freqCost.map;

%% === Initialize uvBias ===
NNF_H.uvBias.map    = im2single(zeros(imgH_H, imgW_H, 3));
NNF_H.uvBias.data = uvMat_from_uvMap(NNF_L.uvBias.map, uvPixL);
NNF_H.uvBias.map = update_uvMap(NNF_H.uvBias.map, NNF_H.uvBias.data, NNF_H.uvPix, true(1,NNF_H.uvPix.numUvPix));

%% === Initialize uvCost ===
NNF_H.uvCost.map  = im2single(zeros(imgH_H, imgW_H));
NNF_H.uvCost.data = uvMat_from_uvMap(NNF_L.uvCost.map, uvPixL);
NNF_H.uvCost.data(:,:) = 0;
NNF_H.uvCost.map = update_uvMap(NNF_H.uvCost.map, NNF_H.uvCost.data, NNF_H.uvPix, true(1,NNF_H.uvPix.numUvPix));


%% === Initialize pSize weight ===
NNF_H.pSizeWeight.LvlCur = NNF_L.pSizeWeight.LvlCur - 1;    % current level
NNF_H.pSizeWeight.LvlInd = NNF_H.pSizeWeight.LvlCur:min(NNF_H.pSizeWeight.LvlCur+optS.numPSizeLvl-1,optS.numPyrLvl);    % 当前lv下需要处理的哪几个lv
NNF_H.pSizeWeight.data = get_psize_weight(NNF_H.uvPix.dist, pSizehist, NNF_H.pSizeWeight.LvlInd, optS);
NNF_H.pSizeWeight.map  = zeros(NNF_H.imgH, NNF_H.imgW, size(NNF_H.pSizeWeight.LvlInd,2));
NNF_H.pSizeWeight.map = update_uvMap(NNF_H.pSizeWeight.map, NNF_H.pSizeWeight.data, NNF_H.uvPix, true(1,NNF_H.uvPix.numUvPix));

end