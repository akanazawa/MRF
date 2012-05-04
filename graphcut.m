function [E, labels] = graphcut(sites, I, CRF)
%%%%%%%%%%%%%%%%%%%%
% MAP estimate of the CRF cost using the graph-cut algorithm from:
% Multi-label optimization http://vision.csd.uwo.ca/code/ 
% an implementation of alpha-expansion and alpha-beta
% 
% Angjoo Kanazawa 5/3/'12
%%%%%%%%%%%%%%%%%%%%

    L = 2; % binary label for foreground background segmentation
    N = numel(sites);

    h = GCO_Create(N, L); 
    unaryCosts = getUnaryCost(I); %  L by N                                       
    weights = getWeights(sites, CRF, I); % specify the CRF cost here

    % set GCO parameters:
    GCO_SetDataCost(h, unaryCosts);
    % GCO_SetSmoothCost use default Potts model
    GCO_SetNeighbors(h, weights);
    GCO_Expansion(h);
    labels = GCO_GetLabeling(h);
    [E D S L] = GCO_ComputeEnergy(h);
end

%% construct N by N matrix holding unary cost with CRF
function weight = getWeights(sites, CRF, I)
    N = numel(sites);
    weight = zeros(N, N);
    for i = 1:N
        Nx = sites{i}.xNeighbors; 
        Ny = sites{i}.yNeighbors;
        weight(i, Nx) = 4*CRF.K*(CRF.Mx - CRF.Ix(Nx))./CRF.Mx;
        weight(i, Ny) = 4*CRF.K*(CRF.My - CRF.Iy(Ny))./CRF.My;
        % weight(i, Nx) = exp(-CRF.meanDiffSq*(I(Nx)-I(i)).^2)';
        % weight(i, Ny) = exp(-CRF.meanDiffSq*(I(Ny)-I(i)).^2)';
    end
end

