function [U] = pairwise(labels, K, I)
%%%%%%%%%%%%%%%%%%%%
% pairwise.m
% for all labels, compute the pairwise potential 
% B_ij = K if l_i ~= l_j
%        0 else
%
% 
% Angjoo Kanazawa 5/1/'12
%%%%%%%%%%%%%%%%%%%%
[w h] = size(I);
U = 0;
for i = 1:numel(sites)
    N = getNeighbors(i, w, h);
    U = U + K*numel(find(labels(N) ~= labels(x)));
end
