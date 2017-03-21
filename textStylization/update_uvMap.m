function map = update_uvMap(map, data, uvPix, indUpdate)

% UPDATE_UVMAP
%
% Update the map from the data and the indices


uvPixUpdateInd = uvPix.ind(indUpdate);

[imgH, imgW, nCh] = size(map);

for ch = 1: nCh
    map(uvPixUpdateInd + (ch-1)*imgH*imgW) = data(ch,:);
end

end