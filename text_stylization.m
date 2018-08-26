function textEffectFinal = text_stylization(sty, src, srctext, trgtext, optS)

% TEXT_STYLIZATION
%     text effects transfer such that S:S'::T:T'
% 
% Input:
% 	- sty: text effects image
% 	- srctext: source character image
%   - src: the name of the source character (for cache)
% 	- trgtext: target character image
% 	- optS: parameters

% Output:
%	- textEffectFinal: target text effects
% 
% Example:
%   textEffectFinal = text_stylization('flame', 'shu', 'huo', 'imgs/', optS);  
%   
% Disclaimer: 
%   This is a Matlab implementation of the paper:
% 
%   Shuai Yang, Jiaying Liu, Zhouhui Lian and Zongming Guo, 
%   Awesome Typography: Statistics-Based Text Effects Transfer,
%   Accepted by IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2017.
%  
%   This code builds upon the PatchMatch algorithm (https://github.com/jbhuang0604/StructCompletion) 
%   implemented by Jia-Bin Huang.
% 
%   It is provided for educational/research purpose only.
%   If you find the software useful, please consider cite our paper.
%   
%   

% Option parameters
optS = init_opt();

% source text effects S'
srcimg = im2double(sty); 
% source text S
srctext = im2double(rgb2gray(srctext)); 
% target text T
trgtext = im2double(rgb2gray(trgtext));
% target text effects T'
trgimg = zeros(size(trgtext,1), size(trgtext, 2), 3);   
% optinal 
% srcimg = imresize(srcimg, size(trgtext));
% srctext  = imresize(srctext, size(trgtext));

%% === Dist calculation ===

fprintf('- Calculate distance map \n'); 

% get distance map of S and its mean text radius
[srcdist, meandist] = get_dist(srctext, get_sket(srctext, 30));
% get distance map of T and normalize it by the mean text radius of S
trgdist = get_dist(trgtext, get_sket(trgtext, 30), meandist);
% operation on image boundaries
maxdist1 = max([min(trgdist(1,:)),min(trgdist(end,:)),min(trgdist(:,1)), min(trgdist(:,end))]);
maxdist2 = max([min(srcdist(1,:)),min(srcdist(end,:)),min(srcdist(:,1)), min(srcdist(:,end))]);
maxdist = min(maxdist1,maxdist2);
maxdist = 1 * 0.25 + maxdist * 0.75;
srcdist = min(srcdist, maxdist);
trgdist = min(trgdist, maxdist);


%% === Patch size histogram calculation ===

% Construct image pyramid for coarse-to-fine image completion
fprintf('- Construct image pyramid \n'); 

% construct pyramid for target text 
[trgimgPyr, scaleImgPyr] = create_img_pyramid(trgimg, optS);
[trgtextPyr, ~] = create_img_pyramid(trgtext, optS);
[trgdistPyr, ~] = create_img_pyramid(trgdist, optS);
% construct pyramid for source text according to the downsampling scale of the target text 
srctextPyr = create_img_pyramid_by_scale(srctext, optS, scaleImgPyr);
srcimgPyr = create_img_pyramid_by_scale(srcimg, optS, scaleImgPyr);
srcdistPyr = create_img_pyramid_by_scale(srcdist, optS, scaleImgPyr);

fprintf('- Calculate content-aware patch size statistics \n'); 

% get the posterior probability P(l|bin(q))
pSizehist = get_psize_statistics(src, srcimgPyr, scaleImgPyr, srctextPyr, srcdist, optS, true);

%% === Text image stylization ===

% synthesis  
fprintf('- Text effects synthesis \n'); 

% solve objective function 
[trgimgPyr, ~] = synthesis(srcimgPyr, srctextPyr, srcdistPyr, trgimgPyr, trgtextPyr, trgdistPyr, pSizehist, optS);

% return the top level 
textEffectFinal = trgimgPyr{optS.topLevel}; 

end
