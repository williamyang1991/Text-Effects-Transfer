function NNF = init_nnf(trgdistCur, srcdistCur, pSizehist, optS)

% INIT_NNF: Initialize the nearest neighbor field
% Input:
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
[NNF.imgH, NNF.imgW] = size(trgdistCur);
[NNF.imgH2, NNF.imgW2] = size(srcdistCur);
[NNF.validPix, NNF.uvPix] = init_level(trgdistCur, srcdistCur, optS.pRad, optS.lockRandOn, optS.lockRandSeed);

%% === Initialize uvPixN ===
NNF.uvPixN = cell(4,1);
for i = 1: 4
    NNF.uvPixN{i}.sub(1,:) = NNF.uvPix.sub(1,:) - optS.propDir(1,i);
    NNF.uvPixN{i}.sub(2,:) = NNF.uvPix.sub(2,:) - optS.propDir(2,i);
    NNF.uvPixN{i}.ind = sub2ind([NNF.imgH, NNF.imgW], NNF.uvPixN{i}.sub(2,:), NNF.uvPixN{i}.sub(1,:));
    NNF.uvPixN{i}.validInd = NNF.uvPix.mask(NNF.uvPixN{i}.ind);
end

%% === Initialize indMap ===
NNF.pixIndMap = reshape(1:NNF.imgH*NNF.imgW, NNF.imgH, NNF.imgW);
indTrgPatch = prep_target_patch(NNF.pixIndMap, NNF.uvPix.sub, optS);
NNF.trgPatchInd = cat(1, indTrgPatch, ...
    indTrgPatch+NNF.imgH*NNF.imgW, indTrgPatch + 2*NNF.imgH*NNF.imgW);


%% === Initialize uvTform ===
% Initialize unTform.data with random samples
% Distance-based random sampling
if optS.initByDistOn
    randInd = binary_search(NNF.validPix.dist, NNF.uvPix.dist, NNF.validPix.distInd);
% Random sampling
else
    if optS.lockRandOn
        rng(optS.lockRandSeed);
    end
    randInd = randi(NNF.validPix.numValidPix, NNF.uvPix.numUvPix, 1);
end

% matched source patch positions for target patches
uvRandSub = NNF.validPix.sub(:, randInd);
% source patch rotation angles 
randAngle = normrnd(0, 0, 1 , NNF.uvPix.numUvPix);

% init matched source patch (position+rotation)
NNF.uvTform.data = [uvRandSub;randAngle];
NNF.uvTform.map = zeros(NNF.imgH, NNF.imgW, 3);
NNF.uvTform.map = update_uvMap(NNF.uvTform.map, NNF.uvTform.data, NNF.uvPix, true(1,NNF.uvPix.numUvPix));

%% === Initialize source frequency ===
NNF.freq = get_NNF_freq(NNF);
NNF.freqCost.data = zeros(1, NNF.uvPix.numUvPix);
NNF.freqCost.map  = im2single(zeros(NNF.imgH, NNF.imgW));

%% === Initialize uvBias ===
NNF.uvBias.data = zeros(3, NNF.uvPix.numUvPix);
NNF.uvBias.map    = im2single(zeros(NNF.imgH, NNF.imgW, 3));

%% === Initialize uvCost ===
NNF.uvCost.data = zeros(1, NNF.uvPix.numUvPix);
NNF.uvCost.map  = im2single(zeros(NNF.imgH, NNF.imgW));

%% === Initialize update ===
% for early termination
NNF.update.data = false(1, NNF.uvPix.numUvPix);
NNF.update.map  = false(NNF.imgH, NNF.imgW);

%% === Initialize uvPixUpdateSrc ===
NNF.uvPixUpdateSrc.data = zeros(1, NNF.uvPix.numUvPix);
NNF.uvPixUpdateSrc.map  = zeros(NNF.imgH, NNF.imgW);

%% === Initialize pSize weight ===
% in the coarsest level, only one patch scale is used
NNF.pSizeWeight.data = ones(1, NNF.uvPix.numUvPix);
NNF.pSizeWeight.map  = ones(NNF.imgH, NNF.imgW);
NNF.pSizeWeight.LvlInd = optS.numPyrLvl;   
NNF.pSizeWeight.LvlCur = optS.numPyrLvl;   
end