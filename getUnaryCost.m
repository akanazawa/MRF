%% construct L by N matrix holding unary cost 
function [cost, U, labels] = getUnaryCost(sites)
    %%%%%%%%%%
    % Cost definition
    %%%%%%%%%%
    sig = 30;
    mu_b = 128; mu_f1 = 30; mu_f2 = 225;    
    const = 1/2*log(2*pi) + log(sig);
    % alpha_b = @(x) (x - mu_b).^2./(2*sig^2) + const;
    % alpha_f = @(x)-log(exp(-(x - mu_f1).^2./(2*sig^2)) + ...
    %                    exp(-(x - mu_f2).^2./(2*sig^2)) + eps) + ...
    %           const + log(2);    
    %%%%%%%%%%
    alpha_b = (sites - mu_b).^2./(2*sig^2) + const;
    alpha_f = -log(exp(-(sites - mu_f1).^2./(2*sig^2)) + ...
               exp(-(sites - mu_f2).^2./(2*sig^2)) + eps) + ...
               const + log(2);
    cost = [alpha_b, alpha_f]';
    [maxLikelihood argmax] = min(cost);
    labels = argmax - 1;
    U = sum(maxLikelihood);
end
