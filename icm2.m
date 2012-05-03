function [U, labels] = icm2(unary, pairwise, labels)
%%%%%%%%%%%%%%%%%%%%
% icm.m 
% computes the map estimate via Iterated Conditional Modes
% Maximize local conditional probabilities P(f_i | d, f_s/i)
% sequentially until convergence
%
% Angjoo Kanzawa 4/30/'12
%%%%%%%%%%%%%%%%%%%%
    N = size(pairwise, 1);
    % U = getAllEnergy(unary, pairwise, labels);
    fprintf('starting ICM\n');
    numUnchanged = 0;
    while numUnchanged ~= N %diff > 1e-13
        numUnchanged = 0;
        for i = 1:N
            u0 = getEnergy(unary, pairwise, i, labels);
            labels(i) = ~labels(i); %swap
            u1 = getEnergy(unary, pairwise, i, labels);
            if u0 < u1 % no improvement
                labels(i)= ~labels(i); % undo assignment
                numUnchanged = numUnchanged + 1;
            end
        end
        % UNew = getAllEnergy(unary, pairwise, labels);
        % diff = abs(UNew - U);
        % fprintf('\tUold: %g Unew: %g diff: %g\n', U, UNew, diff);
        % U = UNew;
    end
    U = getAllEnergy(unary, pairwise, labels);
end

% get single and pairwise potential for single site
function U = getEnergy(unary, pairwise, ind, labels)
    neigh = find(pairwise(ind, :));
    notSame = find(labels(neigh)~= labels(ind));
    U = unary(labels(ind)+1, ind) + sum(pairwise(ind, neigh(notSame)));
end
