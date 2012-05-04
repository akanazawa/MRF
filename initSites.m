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
    sites{i}.neighbors = getNeighbors(i, [w,h]);
end
end

function [N] = getNeighbors(i, size)
    N = [];
    [x y] = ind2sub(size, i);
    if x + 1 < size(1), N=[N; x+1, y]; end
    if x - 1 > 0, N=[N; x-1, y]; end    
    if y + 1 < size(2), N=[N; x, y+1]; end
    if y - 1 > 0, N=[N; x, y-1]; end    
    N = sub2ind(size, N(:, 2), N(:, 1));
end
