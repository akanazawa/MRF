function sites = initSites(I)
%%%%%%%%%%%%%%%%%%%%
% initSites.m
% makes |I| many cells of struct where each cell has:
%    sites{k}.x = I(i)
%    sites{k}.ind = i
%    site{k}.xNeighbors = [ind1, ind2, ...]
%    site{k}.yNeighbors = [ind1, ind2, ...]
%%%%%%%%%%%%%%%%%%%%

sites = cell(numel(I), 1);
[w, h] = size(I);
for i = 1: numel(I)
    sites{i}.x = I(i);
    sites{i}.ind = i;
    sites{i}.xNeighbors = getX_Neighbors(i, [w,h]);
    sites{i}.yNeighbors = getY_Neighbors(i, [w,h]);    
end

end

function [N] = getNeighbors(i, size)
%%%%%%%%%%%%%%%%%%%%
% getNeighbors.m
% For a given site i=(x,y), get the four neighbors
%
% Angjoo Kanazawa 5/1/'12
%%%%%%%%%%%%%%%%%%%%
    N = [getX_Neighbors(i,size); getY_Neighbors(i,size)];
end

% get neighbors in x-direction
function [N] = getX_Neighbors(i, size)
    [x y] = ind2sub(size, i);
    N = [];
    if x + 1 < size(1), N=[N; x+1, y]; end
    if x - 1 > 0, N=[N; x-1, y]; end    
    N = sub2ind(size, N(:, 2), N(:, 1));
end
% get neighbors in y-direction
function [N] = getY_Neighbors(i, size)
    [x y] = ind2sub(size, i);
    N = [];
    if y + 1 < size(2), N=[N; x, y+1]; end
    if y - 1 > 0, N=[N; x, y-1]; end    
    N = sub2ind(size, N(:, 2), N(:, 1));
end

