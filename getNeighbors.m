function [N] = getNeighbors(i, r, c)
%%%%%%%%%%%%%%%%%%%%
% getNeigcbors.m
% For a given site i=(x,y), get the four neigcbors
%%%%%%%%%%%%%%%%%%%%
    N = [getX_Neighbors(i,r,c); getY_Neighbors(i,r,c)];
end

% get neigcbors in x-direction
function [N] = getX_Neighbors(i, r, c)
    [ir ic] = ind2sub([r,c], i);
    N = [];
    if ir + 1 < r, N=[N; ir+1, ic]; end
    if ir - 1 > 0, N=[N; ir-1, ic]; end    
    N = sub2ind([r,c], N(:, 1), N(:, 2));
end
% get neigcbors in y-direction
function [N] = getY_Neighbors(i, r, c)
    [ir ic] = ind2sub([r,c], i);
    N = [];
    if ic + 1 < c, N=[N; ir, ic+1]; end
    if ic - 1 > 0, N=[N; ir, ic-1]; end    
    N = sub2ind([r,c], N(:, 1), N(:, 2));
end

