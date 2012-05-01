%%%%%%%%%%%%%%%%%%%%
% MRF Problem set 3 CMSC 828 Spring '12
% 
%
% Angjoo Kanazawa 4/30/'12
%%%%%%%%%%%%%%%%%%%%

I = im2double(imread('test.jpg'));
sites = I(:);

%% 1. Data term only
[initLabels, U] = mapUnary(sites);

%% 2. with Pairwise term and ICM
icmLabels = icm(initLabels, sites);

%% 3. with Simulated Annealing
saLabels = sim_annealing(initlabels, sites);

%% 4. CRF
crfLabels = CRF(initlabels, sites); 
