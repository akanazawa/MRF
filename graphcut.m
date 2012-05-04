function [E, labels] = graphcut(unaryCosts, weights)
%%%%%%%%%%%%%%%%%%%%
% MAP estimate of the CRF cost using the graph-cut algorithm from:
% Multi-label optimization http://vision.csd.uwo.ca/code/ 
% an implementation of alpha-expansion and alpha-beta
% 
% Angjoo Kanazawa 5/3/'12
%%%%%%%%%%%%%%%%%%%%
    L = 2; % binary label for foreground background segmentation
    N = size(weights, 1);
    h = GCO_Create(N, L); 
    % set GCO parameters:
    GCO_SetDataCost(h, unaryCosts);
    % GCO_SetSmoothCost use default Potts model
    GCO_SetNeighbors(h, weights);
    GCO_Expansion(h);
    labels = GCO_GetLabeling(h);
    [E D S L] = GCO_ComputeEnergy(h);
end
