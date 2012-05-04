%% construct N by N matrix holding binary cost with/without CRF
function [weight] = getPairwiseCost(I, size, K, CRF)
    N = numel(I);
    weight = sparse(N, N);
    if ~isstruct(CRF)
        for i=1:N
            neigh= getNeighbors(i, size);
            weight(i, neigh) = K;
        end
    else
        for i = 1:N
            [~, Nx, Ny] = getNeighbors(i, size);
            weight(i, Nx) = 4*CRF.K*(CRF.Mx - CRF.Ix(Nx))./CRF.Mx;
            weight(i, Ny) = 4*CRF.K*(CRF.My - CRF.Iy(Ny))./CRF.My;
            % weight(i, Nx) = exp(-CRF.meanDiffSq*(I(Nx)-I(i)).^2)';
            % weight(i, Ny) = exp(-CRF.meanDiffSq*(I(Ny)-I(i)).^2)';
        end
    end
end
function [N, Nx, Ny] = getNeighbors(i, size)
%%%%%%%%%%%%%%%%%%%%
% getNeighbors.m
% For a given site i=(x,y), get the four neighbors
%%%%%%%%%%%%%%%%%%%%
    r= size(1); c = size(2);
    Nx = getX_Neighbors(i, r,c);            
    Ny = getY_Neighbors(i, r,c);                
    N = [Nx; Ny];
end
% get neigcbors in y-direction
function [N] = getY_Neighbors(i, r, c)
    [ir ic] = ind2sub([r,c], i);
    N = [];
    if ir + 1 < r, N=[N; ir+1, ic]; end
    if ir - 1 > 0, N=[N; ir-1, ic]; end    
    N = sub2ind([r,c], N(:, 1), N(:, 2));
end
% get neigcbors in x-direction
function [N] = getX_Neighbors(i, r, c)
    [ir ic] = ind2sub([r,c], i);
    N = [];
    if ic + 1 < c, N=[N; ir, ic+1]; end
    if ic - 1 > 0, N=[N; ir, ic-1]; end    
    N = sub2ind([r,c], N(:, 1), N(:, 2));
end