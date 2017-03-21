function ind = binary_search(L, T, index)

% BINARY_SEARCH: find the nearest element in L for element in T
% Input:
%	- L
%	- T
%	- index: sorting index for T 
% Output:
%   - ind

    ind = zeros(size(T));
    if nargin == 2
        [L, index] = sort(L);
    end
    num = size(L, 2);
    for i = 1:size(T,2)
        target = T(i);
        if target >= L(end)
            ind(i) = num;
            continue;
        elseif target <= L(1)
            ind(i) = 1;
            continue;
        end
        left = 1;
        right = num;
        while right > left
            mid = floor((left + right)/2);
            if L(mid) < target
                left = mid;
            elseif L(mid) > target
                right = mid;
            else
                left = mid;
                break;
            end
            if left == right - 1 && L(left) < target && L(right) > target
                break;
            end
        end
        if abs(L(left)-target) < abs(L(right)-target)
            ind(i) = left;
        else
            ind(i) = right;
        end       
    end
    ind = min(max(1, ind), num);
    ind = index(ind);
end