%% construct N by N matrix holding binary cost with/without CRF
function [weight] = getPairwiseCost(I, neighbors, CRF)
    [r c] = size(I);
    N = numel(I);
    weight = sparse(N, N);
    if ~CRF.isCRF
        for i=1:N
            N = neighbors{i};
            weight(i, N) = CRF.K;
        end
    else
        for i = 1:N
            [Nx, Ny] = getXYNeighbors(i, r, c);
            weight(i, Nx) = 4*CRF.K*(CRF.Mx - CRF.Ix(Nx))./CRF.Mx;
            weight(i, Ny) = 4*CRF.K*(CRF.My - CRF.Iy(Ny))./CRF.My;
            % weight(i, Nx) = exp(-CRF.meanDiffSq*(I(Nx)-I(i)).^2)';
            % weight(i, Ny) = exp(-CRF.meanDiffSq*(I(Ny)-I(i)).^2)';
        end
    end
end
function [Nx, Ny] = getXYNeighbors(i, r, c)
    [ir ic] = ind2sub([r, c], i);
    Nx = []; Ny = [];
    if ir < r, Nx=[Nx; i+1]; end
    if ir > 1, Nx=[Nx; i-1]; end  
    if ic < c, Ny=[Ny; i+r]; end
    if ic > 1, Ny=[Ny; i-r]; end    
end
