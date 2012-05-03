function [U, labels] = mapUnary(sites)
%%%%%%%%%%%%%%%%%%%%
% MAP estimate with unary only
%%%%%%%%%%%%%%%%%%%%

sig = 30;
mu_b = 128; mu_f1 = 30; mu_f2 = 225;

const = 1/2*log(2*pi) + log(sig);
probB = 1/(sqrt(2*pi)*sig)*exp(-(sites - mu_b).^2./(2*sig^2));

alpha_b = (sites - mu_b).^2./(2*sig^2) + const;
% gam1 = exp(-(sites-mu_f1).^2./(2*sig^2))*(sqrt(2*pi)*sig)^(-1);
% gam2 = exp(-(sites-mu_f2).^2./(2*sig^2))*(sqrt(2*pi)*sig)^(-1);
% alpha_f = -gam1 + (sites - mu_f1).^2./(2*sig^2) + ...
%           -gam2 + (sites - mu_f2).^2./(2*sig^2) + ...
%           const + const;

alpha_f = -log(exp(-(sites - mu_f1).^2./(2*sig^2)) + ...
               exp(-(sites - mu_f2).^2./(2*sig^2)) + eps) + ...
               const + log(2);

probF = 1/(2*sqrt(2*pi)*sig)*(exp(-(sites - mu_f1).^2./(2*sig^2)) + exp(-(sites - mu_f2).^2./(2*sig^2)));

[maxLikelihood argmax] = min([alpha_b, alpha_f]');
%[maxLikelihood argmax] = max([probB, probF]');
labels = argmax - 1;
U = sum(maxLikelihood);

end
