function [trgimgPyr, NNF] = one_pass(trgimgPyr, trgtextPyr, trgdistCur, srcimgPyr, srctextPyr, srcdistCur, NNF, wDistPatch, numIterLvl, iLvl, optS, lockAngleFlag)

% ONE_PASS
% 
% Modified PatchMatch algorithm
%
% Update the nearest neighbor field using the PatchMatch algorithm
%
% Input:
%   - trgimgPyr, trgtextPyr, trgdistCur: source of P(p), P'(p), dist(p)
%   - srcimgPyr, srctextPyr, srcdistCur: source of Q(q), Q'(q), dist(q) 
%   - wDistPatch: voting weight for patches
%   - NNF    
%   - wDistPatch
%   - optS
%   - iLvl, numIterLvl: current level, current iteration
%   - lockAngleFlag: whether to use patch rotation
% Output:
%   - NNF
%   - trgimgPyr

img = trgimgPyr{iLvl};
trgtextCur = trgtextPyr{iLvl};
srcimgCur = srcimgPyr{iLvl};
srctextCur = srctextPyr{iLvl};

% for visualization
vissrctextCur =  srctextCur(:,:,ones(3,1));
vistrgtextCur = zeros(size(img));

for iter = 1 : numIterLvl
    
    % === Compute the patch matching cost at the current level ===

    % Prepare target patches  P(p), P'(p), dist(p)
    trgPatchPyr = prep_target_patchPyr(trgimgPyr, NNF.uvPix.sub, NNF.pSizeWeight.LvlInd, iLvl, optS);%prep_target_patch(img, NNF.uvPix.sub,  optS);
    trgTextPatchPyr = prep_target_patchPyr(trgtextPyr, NNF.uvPix.sub, NNF.pSizeWeight.LvlInd, iLvl, optS);
    trgDistPatch = prep_target_patch(trgdistCur, NNF.uvPix.sub, optS);
    for i = 1:size(trgTextPatchPyr, 1)
        trgTextPatchPyr{i} = reshape(trgTextPatchPyr{i}, optS.pNumPix, 1, NNF.uvPix.numUvPix);
    end    
    trgDistPatch = reshape(trgDistPatch, optS.pNumPix, 1, NNF.uvPix.numUvPix);
    
    % Prepare source patches Q(q), Q'(q), dist(q) based on NNF
    srcPatchPyr = prep_source_patchPyr(srcimgPyr, NNF.uvTform.data, NNF.pSizeWeight.LvlInd, iLvl, optS);
    srcTextPatchPyr = prep_source_patchPyr(srctextPyr, NNF.uvTform.data, NNF.pSizeWeight.LvlInd, iLvl, optS);
    srcDistPatch = prep_source_patch(srcdistCur, NNF.uvTform.data, optS);   
    for i = 1:size(srcTextPatchPyr, 1)
        srcTextPatchPyr{i} = reshape(srcTextPatchPyr{i}, optS.pNumPix, 1, NNF.uvPix.numUvPix);
    end 
    srcDistPatch = reshape(srcDistPatch, optS.pNumPix, 1, NNF.uvPix.numUvPix);
   
    % Compute patch matching cost    
    srcInd = sub2ind(size(NNF.validPix.mask), round(NNF.uvTform.data(2,:)), round(NNF.uvTform.data(1,:)));
    [uvCostcur, NNF.uvBias.data] = patch_cost(trgPatchPyr, srcPatchPyr, trgTextPatchPyr, srcTextPatchPyr, trgDistPatch, srcDistPatch, wDistPatch, ...
        NNF.freq.map, srcInd, NNF.pSizeWeight.data, optS, iLvl, iter, numIterLvl);
    NNF.uvCost.data = sum(uvCostcur, 1);
    NNF.uvCost.map = update_uvMap(NNF.uvCost.map, NNF.uvCost.data, NNF.uvPix, true(1, NNF.uvPix.numUvPix));
   
    % Initialize update index map (for early termination)
    NNF.update.data = false(1, NNF.uvPix.numUvPix);
    NNF.update.map  = false(NNF.imgH, NNF.imgW);
    
    % === Update the NNF using the PatchMatch algorithm ===
    % find better matches in four-neighborhood 
    % randomly find better matches within a shrunken window
    [NNF, nUpdate]= update_NNF(trgPatchPyr, srcimgPyr, trgTextPatchPyr, srctextPyr, trgDistPatch, srcdistCur, wDistPatch, NNF, iLvl, optS, iter, numIterLvl, lockAngleFlag);
    avgPatchCost = mean(NNF.uvCost.data, 2);
    
    % === Update the image ===
    img = voting_update(wDistPatch, img, srcimgCur, NNF, optS, optS.useBiasCorrection);
    if(optS.visOn) 
        vistrgtextCur = voting_update(wDistPatch, trgtextCur, srctextCur, NNF, optS);
    end
    trgimgPyr{iLvl} = img;
    for i = optS.numPyrLvl:-1:iLvl+1
        trgimgPyr{i} = imresize(trgimgPyr{iLvl}, [size(trgimgPyr{i},1), size(trgimgPyr{i},2)], optS.resampleKernel);
    end
    
    % === Visualizing the progress ===
    if(optS.visOn) 
        NNFVis = vis_nnf(NNF);
        figure(1024);
        subplot(2,3,1), imshow(img);
        title('Current Image', 'fontsize', 8);
        subplot(2,3,2), imshow(NNFVis.uvTfomMapVis);
        title('Nearest neighbor field', 'fontsize', 8);
        subplot(2,3,3), imshow(vissrctextCur);
        title('Source Text', 'fontsize', 8);      
        subplot(2,3,4), imshow(NNFVis.uvPixUpdateSrcMap);
        title('Update Map', 'fontsize', 8); colormap jet
        subplot(2,3,5), imshow(min(1, NNFVis.uvCostMapVis));
        title('Patch matching cost', 'fontsize', 8); colormap jet
        subplot(2,3,6), imshow(vistrgtextCur);
        title('Target Text', 'fontsize', 8);           
    end
    fprintf('    %3d\t%12d\t%12d\t%12d\t%14f\n', iter, nUpdate(1), nUpdate(2), nUpdate(3), avgPatchCost);
end

% === Visualizing the progress ===
if(~optS.visOn) 
    NNFVis = vis_nnf(NNF);
    figure(1024);
    subplot(2,3,1), imshow(img);
    title('Current Image', 'fontsize', 8);
    subplot(2,3,2), imshow(NNFVis.uvTfomMapVis);
    title('Nearest neighbor field', 'fontsize', 8);
    subplot(2,3,3), imagesc(NNFVis.angleMapVis);axis image
    title('Angle map', 'fontsize', 8);     
    subplot(2,3,4), imshow(NNFVis.uvPixUpdateSrcMap);
    title('Update Map', 'fontsize', 8); 
    subplot(2,3,5), imagesc(min(1, NNFVis.uvCostMapVis));axis image
    title('Patch matching cost', 'fontsize', 8); 
    subplot(2,3,6), imagesc(min(10, NNF.freqCost.map));axis image
    title('Repetitiveness cost', 'fontsize', 8);                    
end

end