function [N] = getNeighbors(i, w, h)
%%%%%%%%%%%%%%%%%%%%
% getNeighbors.m
% For a given site i=(x,y), get the four neighbors
%
% Angjoo Kanazawa 5/1/'12
%%%%%%%%%%%%%%%%%%%%
    N = [getX_Neighbors(i,w,h); getY_Neighbors(i,w,h)];
end

% get neighbors in x-direction
function [N] = getX_Neighbors(i, w, h)
    [x y] = ind2sub([w,h], i);
    N = [];
    if x + 1 < w, N=[N; x+1, y]; end
    if x - 1 > 0, N=[N; x-1, y]; end    
    N = sub2ind([w,h], N(:, 2), N(:, 1));
end
% get neighbors in y-direction
function [N] = getY_Neighbors(i, w, h)
    [x y] = ind2sub([w,h], i);
    N = [];
    if y + 1 < h, N=[N; x, y+1]; end
    if y - 1 > 0, N=[N; x, y-1]; end    
    N = sub2ind([w,h], N(:, 2), N(:, 1));
end

