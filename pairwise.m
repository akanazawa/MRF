function [U] = pairwise(site, labels, K, CRF)
%%%%%%%%%%%%%%%%%%%%
% pairwise.m
% for all labels, compute the pairwise potential 
% B_ij = K if l_i ~= l_j
%        0 else
%
% For CRF, send a struct as option
% If crf is set, 4*K*(M - d/dx(i,j))./M, where M is
% the max derivative in each direction.
% 
% Angjoo Kanazawa 5/1/'12
%%%%%%%%%%%%%%%%%%%%
    Nx = site.xNeighbors; 
    Ny = site.yNeighbors;
    if ~isstruct(CRF)
        U = K*numel(find(labels([Nx ; Ny]) ~= labels(site.ind)));
    else
        indx = find(labels(Nx) ~= labels(site.ind));
        indy = find(labels(Ny) ~= labels(site.ind));
        U = 0;
        if ~isempty(indx); U = 4*K*(CRF.Mx - CRF.Ix(Nx(indx)))./CRF.Mx; end;
        if ~isempty(indy); U = U + 4*K*(CRF.My - CRF.Iy(Ny(indy)))./CRF.My; end;
        U = sum(U);
    end
    
end
