function U = getEnergy(sites, labels, K)
%%%%%%%%%%%%%%%%%%%%
% getEnergy.m
% gets energy of labels
%
% Angjoo Kanazawa 5/1/'12
%%%%%%%%%%%%%%%%%%%%

U = unary(sites) + pairwise(labels, K, I);




