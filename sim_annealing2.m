function [U, labels] = sim_annealing2(unary, pairwise, labels)
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
        for i = 1:numItr
            index = randi(N, 1);
            labelprime = labels;
            labelprime(index) = ~labels(index);
            pf = getProb(unary, pairwise, index, labels, T);
            pfprime = getProb(unary, pairwise, index, labelprime, T);
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

% get p(f_i), the probability of labeling site i with f_i
function p = getProb(unary, pairwise, ind, labels, T)
    neigh = find(pairwise(ind, :));
    notSame = find(labels(neigh)~= labels(ind));
    U = unary(labels(ind)+1, ind) + sum(pairwise(ind, neigh(notSame)));
    p = exp(-U/T);
end

