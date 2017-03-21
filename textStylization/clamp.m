function v = clamp(v, vLowerB, vUpperB)

% CLAMP: clamp value v

v(v<vLowerB) = vLowerB;
v(v>vUpperB) = vUpperB;

end