function pSizeWeight = get_psize_weight(dist, pSizehist, LvlInd, optS)

% GET_PSIZE_WEIGHT
%
% get weight for adaptive scale-aware patch distance at the current level
%
% Input: 
%   - dist
%   - pSizehist
%   - LvlInd: scale from l~min(10,l+L-1)
%   - optS
% Output:
%   - pSizeWeight
%  

% get bin(p) for each p
maxDist = max(dist);
HistNum = optS.histNum;
dist = round((HistNum-1) * dist / maxDist)+1;
dist = max(dist, 1);
dist = min(dist, HistNum);

% posterior probability of P(1|bin(p))~P(10|bin(p)) for each p
pSizehist = pSizehist(dist,:);

% posterior probability of P(l|bin(p))~P(min(10,l+L-1)|bin(p)) for each p
% at level l, the scale used are from l to min(10,l+L-1), rather from 1 to 10
% recompute the weight to make the sum of the posterior probability equals to one
pSizeWeight = pSizehist(:,LvlInd);
weightSum = sum(pSizeWeight,2);
% if current scales do not cover the main part of the probability
Ind = weightSum>=0.3;
if LvlInd(end) == optS.numPyrLvl
    % large scale take the main part    
    pSizeWeight(~Ind,1) = 1;    
elseif LvlInd(1) == 1
    % small scale take the main part    
    pSizeWeight(~Ind,end) = 1;
else
    % compare large scale and small scale
    Ind2 = sum(pSizehist(:,1:LvlInd(1)-1),2) > sum(pSizehist(:,LvlInd(1)+1:end),2);
    pSizeWeight((~Ind)&Ind2,1) = 1;
    pSizeWeight((~Ind)&(~Ind2),end) = 1;
end
% recompute the weight
pSizeWeight = bsxfun(@times, pSizeWeight, 1./sum(pSizeWeight,2));
pSizeWeight = pSizeWeight';

if sum(abs(sum(pSizeWeight,1)-1)>0.001) ~= 0
    fprintf('ERROR: PSIZEWEIGHT NOT SUM TO ONE!!!!\n');
end