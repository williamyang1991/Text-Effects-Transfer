function [distMap, meandist] = get_dist(imgtext, imgsket, meandist)

% GET_DIST
%	get the normalized distance of the text image
% 
% Input:
% 	- imgtext: text image
% 	- imgsket: skeleton
% 	- meandist: (optinal) mean text radius
% Output:
%	- distMap: normalized distance
%	- meandist: mean text radius


I = im2bw(imgtext(:,:,1));
[h, w] = size(I);

% find text contour
edge=bwperim(I);

% distance to the text contour
[distMap1, imap] = bwdist(edge, 'euclidean');
% distance to the text skeleton
[distMap2] = bwdist(imgsket, 'euclidean');
imap = imap(:);
distMap2 = distMap2(:);
% r(q©Ø)
distMap3 = reshape(distMap2(imap), h, w);

% r(q)
py = sort(distMap2(edge));
% rank(q)
px = (round(size(py)*0.2):round(size(py)*0.8))';
% linear regression
p=polyfit(px,py(px),1);
% eliminate outliers (Eq. (5) in the paper)
distMap3 = max(distMap3, size(py,1)*0.2*p(1)+p(2));

% mean text radius
if nargin == 2
    meandist = (size(py,1)*0.5*p(1)+p(2));
end

% normalize distance (Eq. (6) in the paper)
distMap1(I) = -distMap1(I);
distMap = (distMap1 + distMap3) ./ distMap3;
distMap = max(distMap, 0);  % clamp
distMap1 = distMap1 / meandist + 1;
distMap(~I) = distMap1(~I);


