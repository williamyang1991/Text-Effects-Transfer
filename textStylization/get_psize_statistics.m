function pSizehist = get_psize_statistics(srcName, imgPyr, scalePyr, textPyr, dist, optS, ifvisual)

% GET_PSIZE_STATISTICS
%     get optimal scale of each patch and get the posterior probability P(l|bin(q))
% 
% Input:
% 	- srcName
% 	- imgPyr
% 	- scalePyr
% 	- textPyr
% 	- dist: distance map
% 	- optS
%   - ifvisual: if show result
% Output:
%	- pSizehist: posterior probability P(l|bin(q))

%% dist-patchsize-histogram statistics

if nargin == 6
    ifvisual = false;
end

scaleMapFileName = ['cache/psize-dist-stat/',srcName, '-scalemap-',num2str(optS.pSize),'-',num2str(optS.numPyrLvl),'-',num2str(optS.histNum), '.png'];
histFileName = ['cache/psize-dist-stat/',srcName, '-pSizehist-',num2str(optS.pSize),'-',num2str(optS.numPyrLvl),'-',num2str(optS.histNum), '.png'];

% read the cache if have been stored
if(exist(fullfile(histFileName), 'file'))
    pSizehist = im2double(imread(histFileName));   
    return;
end


HistNum = optS.histNum;

% get patch centers q
uvPix = get_pixel_ind(textPyr{1}, optS.pRad);

% optimal scales for q
scaleIndex = ones(1, uvPix.numUvPix);
% patches that satisfy the filter criterion
updateIndex = true(1, uvPix.numUvPix);

fprintf('--- Processing');

% no cache
if(~exist(fullfile(histFileName), 'file'))
    % from scale L to 1
    for i = optS.numPyrLvl:-1:2
        fprintf('.');
        
        % prepare source patches Q(q)
        srcuvPix = get_pixel_ind(imgPyr{i}, optS.pRad);
        srcpatches = prep_target_patch(imgPyr{i}, srcuvPix.sub,  optS);
        srctextpatches = prep_target_patch(textPyr{i}, srcuvPix.sub,  optS);
        srcpatches = reshape(srcpatches, size(srcpatches,1)*size(srcpatches,2),size(srcpatches,3));    
        srcpatches = [srcpatches; srctextpatches];
        
        % prepare target patches Q(q-hat) 
        patches = prep_target_patch(imgPyr{i}, uvPix.sub(:,updateIndex)*scalePyr{i}.imgScale,  optS);
        textpatches = prep_target_patch(textPyr{i}, uvPix.sub(:,updateIndex)*scalePyr{i}.imgScale,  optS);    
        patches = reshape(patches, size(patches,1)*size(patches,2),size(patches,3));      
        patches = [patches; textpatches];
        sigma = std(patches).^(0.5);  
        
        % FLANN search q-hat for q (Eq. (3) in the paper) 
        build_params.target_precision = 1;
        build_params.build_weight = 0.5;
        build_params.memory_weight = 0;
        [index, parameters] = flann_build_index(srcpatches, build_params); 
        [idx, dists] = flann_search(index, patches, 6, parameters);
        
        % FLANN can be replaced by knnsearch
        % However, knnsearch is much slower
        % [idx, dists] = knnsearch(srcpatches', patches', 'K', 6);
        % dists = dists';
        
        % filter criterion (Eq. (4) in the paper) 
        stopIndex = 0.5*sigma + dists(6,:) <= 0.3;
        udpateIndexCur = find(updateIndex==1);
        updateIndex(udpateIndexCur(stopIndex)) = false;
        scaleIndex(udpateIndexCur(stopIndex)) = i;
    end
    fprintf('\n');
    scalemap = zeros(size(textPyr{1}));
    scalemap = update_uvMap(scalemap, scaleIndex, uvPix, true(1, uvPix.numUvPix));
    imwrite(uint8(scalemap), scaleMapFileName);
end

% visualzation (Fig. 3 in the paper)
map = zeros(size(textPyr{1}));
colormap = jet(11);
if ifvisual
    figure, imagesc(scalemap);axis image;
    figure, imshow(imgPyr{1}*0.5);hold on;
    maxpsize = round(optS.pSize/scalePyr{optS.numPyrLvl}.imgScale/2);    
    for i = maxpsize+1:size(textPyr{1},1)-maxpsize-1;
        for j = maxpsize+1:size(textPyr{1},2)-maxpsize-1;
            if scalemap(i, j) == 0 || map(i, j) == 1
                continue;
            end                       
            psize = round(optS.pSize/scalePyr{scalemap(i, j)}.imgScale);
            halfpsize = floor(psize / 2);
            if sum(map(i-halfpsize:i+halfpsize,j-halfpsize:j+halfpsize)) <= halfpsize                   
                plot(j, i, 'LineWidth', 1, 'Marker', 's', 'MarkerSize', psize, 'MarkerEdgeColor', colormap(1+scalemap(i, j),:));        
                if (scalemap(i, j) == 1)
                    map(i-halfpsize*2:i+halfpsize*2,j-halfpsize*2:j+halfpsize*2) = 1;
                elseif (scalemap(i, j) == 10)
                    map(i-halfpsize:i+halfpsize,j-halfpsize:j+halfpsize) = 1;    
                else
                    map(i-halfpsize/2:i+halfpsize/2,j-halfpsize/2:j+halfpsize/2) = 1;
                end
                
            end
        end
    end
end


% get posterior probability P(l|bin(q))
pSizehist = zeros(HistNum, optS.numPyrLvl);

distdata = dist(1+optS.pRad:end-optS.pRad, 1+optS.pRad:end-optS.pRad);
scaledata =  scalemap(1+optS.pRad:end-optS.pRad, 1+optS.pRad:end-optS.pRad);
distdata = round(HistNum * distdata(:) / max(max(distdata)));
distdata = [distdata-1,distdata,distdata+1];
distdata = min(max(distdata, 1), HistNum);
scaledata = scaledata(:);

for i = 1:size(distdata,1)
    pSizehist(distdata(i,1),scaledata(i)) = pSizehist(distdata(i,1),scaledata(i)) + 1;
    pSizehist(distdata(i,2),scaledata(i)) = pSizehist(distdata(i,2),scaledata(i)) + 2;
    pSizehist(distdata(i,3),scaledata(i)) = pSizehist(distdata(i,3),scaledata(i)) + 1;
end

pSizehist = pSizehist ./  repmat(sum((pSizehist'))',1,optS.numPyrLvl);

if ifvisual
    figure, mesh(pSizehist);
    figure, bar3(pSizehist);
end

% cache the result
imwrite(pSizehist, histFileName);