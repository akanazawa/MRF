function [U, labels] = sim_annealing(sites, labels, K, beta, CRF)
%%%%%%%%%%%%%%%%%%%%
% sim_annealing.m
% Minimize the energy by simulated annealing
%
%%%%%%%%%%%%%%%%%%%%    
    temp = [64 32 16 8 4 2 1 .5 .25 .125 .0625 .03125 eps];
    numItr = 81920;
    N = numel(sites);
    randInd = randperm(numel(sites));
    U = getAllEnergy(sites, labels, K, beta, CRF);
    fprintf('starting simulated annealing\n');
    for T = temp
        for i = 1:numItr
            index = randi(N, 1);
            labelprime = labels;
            labelprime(index) = ~labels(index);
            pf = getProb(sites{index}, labels, K, beta, CRF, T);
            pfprime = getProb(sites{index}, labelprime, K, beta, CRF, T);
            % P(f) indicates probability of labeling f
            if rand(1,1) < pfprime/pf
                labels(index) = ~labels(index);
            end
        end
        U2 = getAllEnergy(sites, labels, K, beta, CRF);
        fprintf('\tat T=%d Uold: %g Unew: %g diff: %g\n', ...
                T, U, U2, abs(U-U2));
        U = U2;
    end    
    U = getAllEnergy(sites, labels, K, beta, CRF);
end

% get single and pairwise potential for single site
function p = getProb(site, labels, K, beta, CRF, T)
    U = unary(site, labels(site.ind)) + ...
        beta.*pairwise(site, labels, K, CRF);
    p = exp(-U/T);
end

% get single and pairwise potential for all sites
function U = getAllEnergy(sites, labels, K, beta, CRF)
    U = 0;
    for i = 1:numel(sites)
        U = U + unary(sites{i}, labels(i)) + ...
            beta.*pairwise(sites{i}, labels, K, CRF);
    end
end

