function [freq, freqCost]= get_NNF_freq(NNF)

% GET_NNF_FREQ:
%
% get |¦µ(q)| and the cost of Psycho-Visual Term
%
% Input:
%   - NNF
% Output:
%   - freq: 
%       - freq.data:    1 x numValidPix   
%       - freq.map:     H2 x W2 x 1
%   - freqCost: 
%       - freqCost.data:1 x numUvPix   
%       - freqCost.map: H x W x 1
% 

uvRandSub = round(NNF.uvTform.data(1:2,:));

ind = sub2ind(size(NNF.validPix.mask), uvRandSub(2,:), uvRandSub(1,:));
srcFreqData = zeros(1, NNF.imgH2*NNF.imgW2);
for i = 1:size(ind,2)
    srcFreqData(ind(i)) = srcFreqData(ind(i)) + 1;
end

% from the perspective of q
freq.data = srcFreqData(NNF.validPix.ind);
freq.map = zeros(NNF.imgH2, NNF.imgW2);
freq.map = update_uvMap(freq.map, freq.data, NNF.validPix, true(1, NNF.validPix.numValidPix));

% from the perspective of p
freqCost.data = srcFreqData(ind);
freqCost.map = zeros(NNF.imgH, NNF.imgW);
freqCost.map = update_uvMap(freqCost.map, freqCost.data, NNF.uvPix, true(1, NNF.uvPix.numUvPix));
