function  [NNF, nUpdateTotal] = random_search(trgPatchPyr, srcimgPyr, trgTextPatchPyr, srcTextPyr, trgDistPatch, srcDist, wDistPatch,  NNF, optS, iLvl, iter, numIterLvl, lockAngleFlag)

% SC_RANDOM_SEARCH: update the nearest neighbor using random sampling
% 给每个块，它分配的变换矩阵T，
% 随机取一个偏移量给它，然后计算新的变换矩阵T'，如果T'比T的cost更低，则更新
% 此随机采样过程共重复optS.numRegSample=5次
% 其中偏移量的取值范围会随着迭代次数呈指数减小

numPSizeLvl = size(trgPatchPyr, 1);

[imgH, imgW, nCh] = size(NNF.uvTform.map);
imgSize = max(imgH, imgW);

uvPix = NNF.uvPix;
numUvPix = size(uvPix.sub, 2);

searchRad = max(imgH, imgW)/2;
if lockAngleFlag
    angleRad = 0;
else
    angleRad = optS.angleRad;
end

nUpdateTotal = 0;

uvPixActiveInd = true(1, numUvPix);
iteration = 1;
while(iteration <= optS.numRandSample)
    % while(searchRad > 1)
    iteration = iteration + 1;
    % Reduce search radius by half
    % 偏移量的范围
    searchRad = searchRad/2;
    angleRad = angleRad/2;
    %     if(searchRad < 1)
    %         break;
    %     end
    % 提取出未知块的变换矩阵，9*N
    uvTformCandCur = uvMat_from_uvMap(NNF.uvTform.map, uvPix);
    
    % Draw random samples
    % 随机选取偏移量获得变换矩阵
    srcPos =zeros(size(uvTformCandCur));
    if optS.lockRandOn
        rng(optS.lockRandSeed);
    end
    srcPos(1:2,:) = uvTformCandCur(1:2,:) + 2*searchRad*(rand(2, numUvPix) - 0.5);
    if optS.lockRandOn
        rng(optS.lockRandSeed);
    end
    srcPos(3,:) = uvTformCandCur(3,:) + 2*angleRad*(rand(1, numUvPix)-0.5);
    %srcPos = round(srcPos);
    srcPos(1,:) = clamp(srcPos(1,:), optS.pRad+1, imgW - optS.pRad);
    srcPos(2,:) = clamp(srcPos(2,:), optS.pRad+1, imgH - optS.pRad);
    srcPos(3,:) = clamp(srcPos(3,:), -optS.angleRad, optS.angleRad);
      
    % Draw plane ID candidate
    %uvPlaneIDCand = sc_draw_plane_id(NNF.uvPlaneID.planeProbAcc);
    
    % Estimate the domain transformation
    % 计算新的变换矩阵T'
    uvTformCand = src_domain_tform(srcPos, NNF.uvPix.sub, 1);
    
    %A = isnan(uvTformCand);
    %if(sum(A(:))~=0)
        %error('uvTformCand contains NaN');
    %end
    
    % === Reject invalid samples ===
    % Check if the scale of the source patch valid
    % 如果放大缩小太过，则舍弃
    %uvTformScale = sc_scale_tform(uvTformCand);
    %uvValidScaleInd = (uvTformScale > optS.minScale) & (uvTformScale < optS.maxScale);
    % Check if the souce patch is valid
    % 如果不是已知块则舍弃    
    uvValidSrcInd = check_valid_uv(uvTformCand(1:2,:), NNF.validPix.mask);
    % Check if the cost is already low
    % 如果当前块的cost已经很小了，就不用更新了    
    uvValidCostInd = NNF.uvCost.data > optS.rsThres;
   
    % 更新的候选集    
    uvValidInd = uvPixActiveInd & uvValidSrcInd & uvValidCostInd;
    
    uvPixActivePos = find(uvValidInd);
    numActPix = size(uvPixActivePos, 2);
    
    % 如果候选集不为空    
    if(numActPix~=0)
        
        % Update
        % 取得候选集
        trgPatchPyrCur = cell(numPSizeLvl, 1);
        trgTextPatchPyrCur = cell(numPSizeLvl, 1);   
        for i = 1:numPSizeLvl
            trgPatchPyrCur{i} = trgPatchPyr{i}(:,:, uvValidInd);
            trgTextPatchPyrCur{i} = trgTextPatchPyr{i}(:,:, uvValidInd);
            %trgPatchCur      = trgPatch(:,:,uvValidInd);
            %trgTextPatchCur  = trgTextPatch(:,:,uvValidInd);
        end        
        
        trgDistPatchCur  = trgDistPatch(:,:,uvValidInd);
        wDistPatchCur    = wDistPatch(:, uvValidInd);
        pSizeWeightCur   = NNF.pSizeWeight.data(:, uvValidInd);
        uvCostDataCur    = NNF.uvCost.data(:,uvValidInd);
        uvTformCandCur   = uvTformCand(:, uvValidInd);
        %uvPlaneIDCandCur = uvPlaneIDCand(uvValidInd);
        
        uvPixValid.sub = uvPix.sub(:,uvValidInd);
        uvPixValid.ind = uvPix.ind(uvValidInd);
        
        %uvDtBdPixPosCur = NNF.uvDtBdPixPos(:, uvValidInd);
          
        % Grab source patches
        
        srcPatchPyr = prep_source_patchPyr(srcimgPyr, uvTformCandCur, NNF.pSizeWeight.LvlInd, iLvl, optS);
        srcTextPatchPyr = prep_source_patchPyr(srcTextPyr, uvTformCandCur, NNF.pSizeWeight.LvlInd, iLvl, optS);
        srcDistPatch = prep_source_patch(srcDist, uvTformCandCur, optS);
        %subplot(2,1,2), imshow(srcTextPatch);

        for i = 1:numPSizeLvl
            srcTextPatchPyr{i} = reshape(srcTextPatchPyr{i}, optS.pNumPix, 1, size(srcDistPatch,2));
        end 
        srcDistPatch = reshape(srcDistPatch, optS.pNumPix, 1, size(srcDistPatch,2));  
                        
        %srcPatch = prep_source_patch(srcimgCurPyr, uvTformCandCur, optS);
        %srcTextPatch = prep_source_patch(srcTextPyr, uvTformCandCur, optS);
        %srcTextPatch = reshape(srcTextPatch, optS.pNumPix, 1, size(srcTextPatch,2));
        %srcDistPatch = prep_source_patch(srcDist, uvTformCandCur, optS);
        %srcDistPatch = reshape(srcDistPatch, optS.pNumPix, 1, size(srcDistPatch,2));        
        
        
        % 使用T'重新计算匹配块
        srcInd = sub2ind(size(NNF.validPix.mask), round(uvTformCandCur(2,:)), round(uvTformCandCur(1,:)));
        
        [costPatchCandAll, uvBiasCand] = ...
            patch_cost(trgPatchPyrCur, srcPatchPyr, trgTextPatchPyrCur, srcTextPatchPyr, trgDistPatchCur, srcDistPatch, ...
            wDistPatchCur, NNF.freq.map, srcInd, pSizeWeightCur, optS, iLvl, iter, numIterLvl);
        costPatchCand = sum(costPatchCandAll, 1);
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
        updateInd = (costPatchCand < uvCostDataCur);
        nUpdate = sum(updateInd);
        
        % 如果有可以更新的块        
        if(nUpdate~=0)
            
            uvPixActivePos = uvPixActivePos(updateInd);
            
            uvPixActiveInd(uvPixActivePos) = 0;
            
            nUpdateTotal = nUpdateTotal + nUpdate;
            
            % === Update NNF data ===
            % 下面就是各种更新NNF的内容了               
            NNF.uvTform.data(:, uvPixActivePos) = uvTformCandCur(:,updateInd);
            %NNF.uvPlaneID.data(uvPixActivePos)  = uvPlaneIDCandCur(updateInd);
            NNF.uvCost.data(uvPixActivePos)     = costPatchCand(updateInd);
            
            NNF.freqCost.data(uvPixActivePos) = costRetCand(updateInd);          
            
            if(optS.useBiasCorrection)
                NNF.uvBias.data(:,uvPixActivePos)   = uvBiasCand(:,updateInd);
            end
            NNF.update.data(uvPixActivePos) = 2;
            
            %
            NNF.uvPixUpdateSrc.data(uvPixActivePos) = 2;
            
            % === Update NNF map ===
            NNF.uvTform.map = update_uvMap(NNF.uvTform.map, uvTformCandCur(:,updateInd), uvPixValid, updateInd);
            %NNF.uvPlaneID.map = sc_update_uvMap(NNF.uvPlaneID.map, uvPlaneIDCandCur(updateInd), uvPixValid, updateInd);
            NNF.uvCost.map  = update_uvMap(NNF.uvCost.map, costPatchCand(updateInd), uvPixValid, updateInd);
            if(optS.useBiasCorrection)
                NNF.uvBias.map  = update_uvMap(NNF.uvBias.map, uvBiasCand(:,updateInd), uvPixValid, updateInd);
            end
            
            NNF.freqCost.map =  update_uvMap(NNF.freqCost.map, costRetCand(updateInd), uvPixValid, updateInd);   
            
            NNF.update.map  = update_uvMap(NNF.update.map, 1, uvPixValid, updateInd);
            NNF.uvPixUpdateSrc.map  = update_uvMap(NNF.uvPixUpdateSrc.map, 2, uvPixValid, updateInd);
                        
        end      
    end
    
[NNF.freq, freqCost] = get_NNF_freq(NNF);
distweight = max(1, abs(NNF.uvPix.dist));
NNF.uvCost.data = NNF.uvCost.data + (freqCost.data - NNF.freqCost.data)*optS.lambdaRep*optS.repCostRatio./distweight;
NNF.uvCost.map  = update_uvMap(NNF.uvCost.map, NNF.uvCost.data , NNF.uvPix, true(1, NNF.uvPix.numUvPix));
%NNF.uvCost.map = NNF.uvCost.map + (freqCost.map - NNF.freqCost.map)*optS.lambdaRep*optS.repCostRatio./distweight;
NNF.freqCost = freqCost;
    
    
end

end