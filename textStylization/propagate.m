function [NNF, nUpdateTotal] = propagate(trgPatchPyr, srcimgPyr, trgTextPatchPyr, srcTextPyr, trgDistPatch, srcDist, wDistPatch, NNF, optS, indDirection, iLvl, iter, numIterLvl)

% SC_PROPAGATE: update the nearest neighbor field using propagation

% Input:
%   - trgPatch, wDistPatch, img, NNF, modelPlane, optS, indDirection
% Output:
%   - NNF, nUpdateTotal
% 对含有未知像素点的目标块，当前分配的变换矩阵为T
% 查看它邻域分配的变换矩阵Tn分配给它之后（记为T'），cost是否下降
% 如果cost减少了，则更新其变换矩阵为T'


%[imgH, imgW, nCh] = size(srcimgCur);
%imgSize = max(imgH, imgW);
numPSizeLvl = size(trgPatchPyr, 1);


% 已经更新了的块的个数
nUpdateTotal = 0;

% The positions of neighboring pixels
% 邻域的块
uvPixN = NNF.uvPixN{indDirection};
% 邻域的有效块(是未知块，带有偏移量的)
uvPixActiveInd = true(1, NNF.uvPix.numUvPix);
uvPixActiveInd = uvPixActiveInd & uvPixN.validInd;

numUpdatePix = NNF.uvPix.numUvPix;
% 如果上一次更新的像素个数不为0，则继续更新
while(numUpdatePix ~= 0)
    % Prepare uvPix, uvPixNCur
    % 有效块
    uvPix.sub     = NNF.uvPix.sub(:, uvPixActiveInd); uvPix.ind     = NNF.uvPix.ind(:, uvPixActiveInd);
    % 对应的邻域块
    uvPixNCur.sub = uvPixN.sub(:, uvPixActiveInd);    uvPixNCur.ind = uvPixN.ind(:, uvPixActiveInd);
    % 块的最近已知块的距离
    %uvDtBdPixPosCur = NNF.uvDtBdPixPos(:, uvPixActiveInd);
    % 有效块的像素值
    trgPatchPyrCur = cell(numPSizeLvl, 1);
    trgTextPatchPyrCur = cell(numPSizeLvl, 1);
    for i = 1:numPSizeLvl
        trgPatchPyrCur{i} = trgPatchPyr{i}(:,:, uvPixActiveInd);
        trgTextPatchPyrCur{i} = trgTextPatchPyr{i}(:,:, uvPixActiveInd);
        %trgPatchCur   = trgPatch(:,:, uvPixActiveInd);
        %trgTextPatchCur   = trgTextPatch(:,:, uvPixActiveInd);
    end
    trgDistPatchCur   = trgDistPatch(:,:, uvPixActiveInd);
    % 块的权重
    wDistPatchCur = wDistPatch(:,uvPixActiveInd);
    pSizeWeightCur = NNF.pSizeWeight.data(:,uvPixActiveInd);
    % 当前分配的偏移量
    srcPosCur     = NNF.uvTform.data(:, uvPixActiveInd);    
    uvCostCur     = NNF.uvCost.data(:, uvPixActiveInd);
    %uvPlaneIDCur  = NNF.uvPlaneID.map(uvPixNCur.ind);
    %uvPixUpdateSrc= NNF.uvPixUpdateSrc.data(uvPixActiveInd);
    
    % Active pixel positions
    % 有效块的坐标
    uvPixActivePos = find(uvPixActiveInd);
    
    % Get candidate uvTform candidates
    % 有校块的邻域分配的变换矩阵Tn
    % uvTformCand: 2*N
    uvTformCand = uvMat_from_uvMap(NNF.uvTform.map, uvPixNCur);
    
    % Generate candidate transformation by propagation
    % 这些邻域的变换矩阵应用到当前有效块上后的变换矩阵T'
    % 即将邻域的变换矩阵用到有效块上
    % 以一维为例，假设indDirection对应的是1
    % 则当前像素位置为x，当前处理的邻域位置为x-1
    % uvTformCand中x-1处的source patch位于y
    % 那么uvTformCand在x处的source patch位于y+1
    %
    % trg     src
    % x-1 -->  y
    %  x  --> y+1
    %
    uvTformCand = trans_tform(uvTformCand, optS.propDir(:,indDirection));
        
    % Check if the nearest neighbors are valid source patches
    % T'变换后的块不是已知块，则无效
    uvValidSrcInd = check_valid_uv(uvTformCand(1:2,:), NNF.validPix.mask);
    % Check if the nearest neighbors are already the same as the existing one
    % T'如果和T相同，就不用更新了
    diff = abs(uvTformCand - srcPosCur);
    uvValidDistInd = ((diff(1,:) > 1 ) | (diff(2,:) > 1 ) | (diff(3,:) > 0.05));
    
    % Valid pixel indices
    % 有效的T'的位置
    uvValidInd = uvValidSrcInd & uvValidDistInd;
    
    numUvValid = sum(uvValidInd);
    
    if(numUvValid ~= 0)
        % 去掉无效的T'位置上的块，保留有效的块
        for i = 1:numPSizeLvl
            trgPatchPyrCur{i} = trgPatchPyrCur{i}(:,:, uvValidInd);
            trgTextPatchPyrCur{i}= trgTextPatchPyrCur{i}(:,:, uvValidInd);
        end
        %trgPatchCur    = trgPatchCur(:,:, uvValidInd);
        %trgTextPatchCur= trgTextPatchCur(:,:, uvValidInd);
        trgDistPatchCur= trgDistPatchCur(:,:, uvValidInd);
        wDistPatchCur  = wDistPatchCur(:,uvValidInd);
        pSizeWeightCur = pSizeWeightCur(:,uvValidInd);
        uvTformCand    = uvTformCand(:, uvValidInd);        
        uvCostCur      = uvCostCur(uvValidInd);
        
        uvPixUpdatePos = uvPixActivePos(uvValidInd);
        uvPixValid.sub = uvPix.sub(:,uvValidInd);
        uvPixValid.ind = uvPix.ind(uvValidInd);
        %uvPlaneIDCand  = uvPlaneIDCur(uvValidInd);
        
        %uvDtBdPixPosCur = uvDtBdPixPosCur(:, uvValidInd);
        
        % Grab source patches
        % 使用T'重新计算匹配块
        
        srcPatchPyr = prep_source_patchPyr(srcimgPyr, uvTformCand, NNF.pSizeWeight.LvlInd, iLvl, optS);
        srcTextPatchPyr = prep_source_patchPyr(srcTextPyr, uvTformCand, NNF.pSizeWeight.LvlInd, iLvl, optS);
        srcDistPatch = prep_source_patch(srcDist, uvTformCand, optS);
        %subplot(2,1,2), imshow(srcTextPatch);

        for i = 1:numPSizeLvl
            srcTextPatchPyr{i} = reshape(srcTextPatchPyr{i}, optS.pNumPix, 1, size(srcDistPatch,2));
        end 
        srcDistPatch = reshape(srcDistPatch, optS.pNumPix, 1, size(srcDistPatch,2));             
                
        % Compute patch matching cost
        % 计算目标块与匹配块的cost值
        srcInd = sub2ind(size(NNF.validPix.mask), round(uvTformCand(2,:)), round(uvTformCand(1,:)));
        [costPatchCandAll, uvBiasCand] = patch_cost(trgPatchPyrCur, srcPatchPyr, trgTextPatchPyrCur, srcTextPatchPyr, trgDistPatchCur, srcDistPatch, ... 
            wDistPatchCur, NNF.freq.map, srcInd, pSizeWeightCur, optS, iLvl, iter, numIterLvl);
                
        costPatchCand = sum(costPatchCandAll, 1);
        %weight = reshape(trgDistPatchCur(optS.pMidPix,:,:), 1, size(trgDistPatchCur, 3));
        % 记录下当前的重复项cost        
        if optS.lambdaRep ~= 0
            freqData = NNF.freq.map(:)';
            costRetCand = freqData(srcInd);            
            %costRetCand = costPatchCandAll(4,:)/optS.lambdaRep/optS.repCostRatio;
        else
            costRetCand = costPatchCandAll(4,:);
        end
        
        % Check which one to update
        % 找到cost更小的
        updateInd = costPatchCand < uvCostCur;
        
        % cost更小的块的位置
        uvPixUpdatePos = uvPixUpdatePos(updateInd);
        % 这些位置上的T是要更新为T'的，这些需要更新的块的数量
        numUpdatePix = size(uvPixUpdatePos, 2);
    else
        numUpdatePix = 0;
    end
   
    % Update NNF data
    % 如果有可以更新的块
    if(numUpdatePix ~= 0)
        nUpdateTotal = nUpdateTotal + numUpdatePix;
        
        % 下面就是各种更新NNF的内容了
        % === Update NNF data ===
        NNF.uvTform.data(:, uvPixUpdatePos) = uvTformCand(:,updateInd);
        NNF.uvCost.data(uvPixUpdatePos)     = costPatchCand(updateInd);        
        %NNF.uvPlaneID.data(uvPixUpdatePos)  = uvPlaneIDCand(updateInd);
        
        NNF.freqCost.data(uvPixUpdatePos) = costRetCand(updateInd);          
        
        % Apply bias correction
        if(optS.useBiasCorrection)
            NNF.uvBias.data(:,uvPixUpdatePos)   = uvBiasCand(:,updateInd);
        end
        NNF.update.data(:,uvPixUpdatePos)   = 1; % updateInd;
        
        % Label as update by propagation
        NNF.uvPixUpdateSrc.data(uvPixUpdatePos) = 3;
        
        % === Update NNF map ===
        % 将2*N转换到m*n*2
        NNF.uvTform.map    = update_uvMap(NNF.uvTform.map, uvTformCand(:,updateInd), uvPixValid, updateInd);
        NNF.uvCost.map     = update_uvMap(NNF.uvCost.map, costPatchCand(updateInd), uvPixValid, updateInd);
        
        NNF.freqCost.map =  update_uvMap(NNF.freqCost.map, costRetCand(updateInd), uvPixValid, updateInd);    
        
        %NNF.uvPlaneID.map  = sc_update_uvMap(NNF.uvPlaneID.map, uvPlaneIDCand(updateInd), uvPixValid, updateInd);
        if(optS.useBiasCorrection)
            NNF.uvBias.map  = update_uvMap(NNF.uvBias.map, uvBiasCand(:,updateInd), uvPixValid, updateInd);
        end
        NNF.update.map  = update_uvMap(NNF.update.map, 1, uvPixValid, updateInd);
        NNF.uvPixUpdateSrc.map = update_uvMap(NNF.uvPixUpdateSrc.map, 3, uvPixValid, updateInd);                 
        
        % === Update uvPixActiveInd ===
        % 以一维为例，假设indDirection对应的是1
        % 则当前像素位置为x，当前处理的邻域位置为x-1
        % x-1处的source patch位于y
        % 已经决定将x-1处的偏移量propagate给x处，则x处的source patch位于y+1
        % 在下一次的更新中，由于此次只更新了x，并且方向为对应的是1
        % 所以下一次只有可能更新x+1处的像素的source patch
        % 下面的uvPixNextInd就是此类x+1的像素点的索引
        
        % trg     src
        % x-1 -->  y
        %  x  --> y+1  (本轮更新了)
        % x+1          (下轮可能会更新，放入uvPixActiveInd)
        uvPixNextSub = uvPixValid.sub(:,updateInd);
        uvPixNextSub(1,:) = uvPixNextSub(1,:) + optS.propDir(1,indDirection);
        uvPixNextSub(2,:) = uvPixNextSub(2,:) + optS.propDir(2,indDirection);
        uvPixNextInd = sub2ind([NNF.imgH, NNF.imgW], uvPixNextSub(2,:), uvPixNextSub(1,:)); % m*n
        
        updateMap = NNF.uvPix.mask; % m*n
        updateMap(uvPixNextInd) = 0;    % 将x+1标记为0
        uvPixActiveInd = ~updateMap(NNF.uvPix.ind); % 1*N，将属于待处理区域的x+1标记为1
        uvPixActiveInd = uvPixActiveInd & uvPixN.validInd;  % 1*N，去掉包含边界区域（未知区域）的点
                
    end
end

% === Update NNF source frequency ===
[NNF.freq, freqCost] = get_NNF_freq(NNF);
distweight = max(1, abs(NNF.uvPix.dist));
NNF.uvCost.data = NNF.uvCost.data + (freqCost.data - NNF.freqCost.data)*optS.lambdaRep*optS.repCostRatio./distweight;
%NNF.uvCost.map = NNF.uvCost.map + (freqCost.map - NNF.freqCost.map)*optS.lambdaRep*optS.repCostRatio./distweight;
NNF.uvCost.map  = update_uvMap(NNF.uvCost.map, NNF.uvCost.data , NNF.uvPix, true(1, NNF.uvPix.numUvPix));
NNF.freqCost = freqCost;

end