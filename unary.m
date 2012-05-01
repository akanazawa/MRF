function [U] = unary(sites, labels)
%%%%%%%%%%%%%%%%%%%%
% unary.m
% for all sites, compute the unary potential given the label
%
% A foreground (l=1) site comes from a mixture of 2 gaussians
% N(30,30), N(225,30)
% A background (l=0) site comes from N(128,30)
%
% i.e. 
% P(x_i|l_i=0) = exp(-(x_i-mu_b)^2/2*sig)/const = exp(-ah_b)
% P(x_i|l_i=1) = 1/2*exp(-(x_i-mu_f1)^2/2*sig)/const + 
%                1/2*exp(-(x_i-mu_f2)^2/2*sig)/const = exp(-ah_f)
%
% so ah_b = (x-mu_b)^2/2sig^2
% ah_f = log( exp((x-mu_f1)^2/2sig^2) + exp((x-mu_f2)^2/2sig^2))
%% where gam is the weighting of how much each gaussian contributes
%
% Angjoo Kanazawa 5/1/'12
%%%%%%%%%%%%%%%%%%%%
sig = 30;
mu_b = 128, mu_f1 = 30, mu_f2 = 225;

const = 1/2*log(2*pi*sig);
alpha_b = @(x)(x - mu_b)^2./2(sig)^2 + const;
alpha_f = @(x)-log(exp(-(x - mu_f1).^2./2(sig)^2 + exp(-(x - mu_f2).^2./2(sig)^2))...
              + const - log(2);
ind_f = find(label==1);
ind_b = find(label==0);
U = alpha_b(sites(ind_f)) + alpha_f(sites(ind_f));

end

function [U, labels] = mapUnary(sites)
%%%%%%%%%%%%%%%%%%%%
%MAP estimate with unary only
%%%%%%%%%%%%%%%%%%%%
sig = 30;
mu_b = 128, mu_f1 = 30, mu_f2 = 225;

const = 1/2*log(2*pi*sig);
alpha_b = (sites - mu_b)^2./2(sig)^2 + const;
% gam = exp(-(sites-mu_f1).^2./2(sig)^2)./sqrt(2*pi)*sig;
% alpha_f = gam*(sites - mu_f1).^2./2(sig)^2 + (1-gam)*(sites - mu_f2).^2./2(sig)^2;
alpha_f = -log(exp(-(sites - mu_f1).^2./2(sig)^2 + exp(-(sites - mu_f2).^2./2(sig)^2))...
              + const - log(2);

[maxLikelihood argmax] = min[alpha_b; alpha_f];
labels = argmax - 1;
U = sum(maxLikelihood);

end
