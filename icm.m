function [U, labels] = icm(sites, labels, I, K, beta, CRF)
%%%%%%%%%%%%%%%%%%%%
% icm.m 
% computes the map estimate via Iterated Conditional Modes
% Maximize local conditional probabilities P(f_i | d, f_s/i)
% sequentially until convergence
%
% Angjoo Kanzawa 4/30/'12
%%%%%%%%%%%%%%%%%%%%
    [w, h] = size(I);
    N = numel(sites);
    % U = getAllEnergy(sites, labels, K, beta, CRF);
    diff = Inf;
    c = 0;
    fprintf('starting ICM\n');
    numUnchanged = 0;
    while  numUnchanged ~= N %diff > 1e-13
        numUnchanged = 0;
        for i = 1:numel(sites)
            u0 = getEnergy(sites{i}, labels, K, beta, CRF);
            labels(i) = ~labels(i); %swap
            u1 = getEnergy(sites{i}, labels, K, beta, CRF);
            if u0 < u1 % no improvement
                labels(i)= ~labels(i); % undo assignment
                numUnchanged = numUnchanged + 1;
            end
        end
          % UNew = getAllEnergy(sites, labels, K, beta, CRF);
          % diff = abs(UNew - U);
          % fprintf('\titer %d Uold: %g Unew: %g diff: %g\n', c, U, UNew, diff);
          % U = UNew;
          % c = c+1;
    end
    U = getAllEnergy(sites, labels, K, beta, CRF);
end

% get single and pairwise potential for single site
function U = getEnergy(site, labels, K, beta, CRF)
    U = unary(site, labels(site.ind)) + ...
        beta.*pairwise(site, labels, K, CRF);
end

% get single and pairwise potential for all sites
function U = getAllEnergy(sites, labels, K, beta, CRF)
    U = 0;
    for i = 1:numel(sites)
        U = U + unary(sites{i}, labels(i)) + ...
            beta.*pairwise(sites{i}, labels, K, CRF);
    end
end

