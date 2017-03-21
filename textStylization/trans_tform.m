function uvTform = trans_tform(uvTformV, d)

% TRANS_TFORM:
%
% update source patch positions after applying the offsets d to it
%
%

sintheta = sin(pi * uvTformV(3,:));
costheta = cos(pi * uvTformV(3,:));

if(size(d, 2) == 1)
    uvTform(1,:) = costheta*d(1) - sintheta*d(2) + uvTformV(1,:);
    uvTform(2,:) = sintheta*d(1) + costheta*d(2) + uvTformV(2,:);
else
    uvTform(1,:) = costheta.*d(1,:) - sintheta.*d(2,:) + uvTformV(1,:);
    uvTform(2,:) = sintheta.*d(1,:) + costheta.*d(2,:) + uvTformV(2,:);
end

uvTform(3,:) = uvTformV(3,:);

end