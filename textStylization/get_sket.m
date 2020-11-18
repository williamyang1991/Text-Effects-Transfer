function sket = get_sket(I, thresh)

% GET_SKET
%	get the skeleton of the text, with small branches are removed.
% 
% Input:
% 	- I: text image
% 	- thresh: the maximal length of the branches to be removed
% Output:
%	- sket: skeleton of I

% rough skeleton
BW = im2bw(I);
sket=bwmorph(BW,'skel', Inf);

% get the connected strokes
[label,n]=bwlabel(sket);

result = im2bw(zeros(size(I)));

% for each connected stroke, remove its small branches
for i = 1:n
    sket = label == i;
    % find its branchpoints
    temp=bwmorph(sket,'branchpoints');
    [pi,pj] = ind2sub(size(sket), find(temp == 1));
    BW2 = sket;
    curthresh = thresh;
    % for each branchpoint
    for k = 1:size(pi)
        i1 = pi(k);
        j1 = pj(k);
        % remove the branchpoint
        if i1 == 1 || j1 == 1 || i1 == size(sket,1) || j1 == size(sket,2)
            continue;
        end
        sket(i1-1:i1+1, j1-1:j1+1) = 0;
        % remove small branches
        sket2 = bwareaopen(sket, thresh);
        % adjust threshold to ensure this connected stroke wil not be all removed
        while(max(max(sket2))==0)
            if curthresh <= 0
                break;
            end
            curthresh = curthresh - 5;
            sket2 = bwareaopen(sket, curthresh);
        end
        sket = sket2;
         % restore the branchpoint
        sket(i1-1:i1+1, j1-1:j1+1) = BW2(i1-1:i1+1, j1-1:j1+1);
    end
    result = max(result, sket);
end

% remove noises
sket = bwareaopen(result, 9);
