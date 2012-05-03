% get single and pairwise potential for all 
function U = getAllEnergy(unary, pairwise, labels)
    N = size(pairwise, 1);
    U = 0;
    for i = 1:N
        neigh = find(pairwise(i, :));
        notSame = find(labels(neigh) ~= labels(i));
        U = U + unary(labels(i)+1, i) + sum(pairwise(i, neigh(notSame))); 
    end
end
