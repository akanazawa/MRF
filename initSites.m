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
[r, c] = size(I);
for i = 1: numel(I)
    sites{i}.x = I(i);
    sites{i}.neighbors = getNeighbors(i, r, c);
end
end

function [N] = getNeighbors(i, r, c)
    N = [];
    [ir ic] = ind2sub([r, c], i);
    if ir + 1 < r, N=[N; ir+1, ic]; end
    if ir - 1 > 0, N=[N; ir-1, ic]; end    
    if ic + 1 < ic, N=[N; ir, ic+1]; end
    if ic - 1 > 0, N=[N; ir, ic-1]; end    
    N = sub2ind([r, c], N(:, 1), N(:, 2));
end
