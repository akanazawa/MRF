function [U, labels] = sim_annealing2(unary, pairwise, labels, sites)
%%%%%%%%%%%%%%%%%%%%
% sim_annealing.m
% Minimize the energy by simulated annealing
%
%%%%%%%%%%%%%%%%%%%%    
    temp = [64 32 16 8 4 2 1 .5 .25 .125 .0625 .03125 eps];
    numItr = 81920;
    N = size(pairwise, 1);
    % U = getAllEnergy(unary, pairwise, labels);
    fprintf('starting simulated annealing\n');
    for T = temp
        inds = randi(N, numItr, 1);
        for i = 1:numItr
            index = inds(i);
            [pf pfprime] = getProb(unary, pairwise, index, labels, sites,T);
            if rand(1,1) < pfprime/pf
                labels(index) = ~labels(index);
            end
        end
        % U2 = getAllEnergy(unary, pairwise, labels);
        % fprintf('\tat T=%d Uold: %g Unew: %g diff: %g\n', ...
        %         T, U, U2, abs(U-U2));
        % U = U2;
    end    
    U = getAllEnergy(unary, pairwise, labels);
end

% get p(f_i), the probability of labeling site i with f_i and
% p(f_i') that of switched label
function [p, pprime] = getProb(unary, pairwise, ind, labels, sites,T)
    neigh = sites{ind}.neighbors; % find(pairwise(ind, :)) takes a while
    notSame = find(labels(neigh)~= labels(ind));
    notSameChanged = find(labels(neigh) == labels(ind));
    u0 = unary(labels(ind)+1, ind) + ...
         sum(pairwise(ind, neigh(notSame)));
    u1 = unary(~labels(ind)+1, ind) + ...
         sum(pairwise(ind, neigh(notSameChanged)));
    p = exp(-u0/T);
    pprime = exp(-u1/T);    
end

