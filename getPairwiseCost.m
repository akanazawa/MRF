%% construct N by N matrix holding binary cost with/without CRF
function weight = getPairwiseCost(I, size, K, CRF)
    N = numel(I);
    weight = zeros(N, N);
    if ~isstruct(CRF)
        for i=1:N
            neigh = getNeighbors(i, size);
            weight(i, neigh) = K;
        end
    else
        for i = 1:N
            Nx = getX_Neighbors(i, size);            
            Ny = getY_Neighbors(i, size);            
            weight(i, Nx) = 4*CRF.K*(CRF.Mx - CRF.Ix(Nx))./CRF.Mx;
            weight(i, Ny) = 4*CRF.K*(CRF.My - CRF.Iy(Ny))./CRF.My;
            % weight(i, Nx) = exp(-CRF.meanDiffSq*(I(Nx)-I(i)).^2)';
            % weight(i, Ny) = exp(-CRF.meanDiffSq*(I(Ny)-I(i)).^2)';
        end
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

