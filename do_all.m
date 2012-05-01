%%%%%%%%%%%%%%%%%%%%
% MRF Problem set 3 CMSC 828 Spring '12
% 
%
% Angjoo Kanazawa 4/30/'12
%%%%%%%%%%%%%%%%%%%%

I = im2double(imread('test.jpg'));
sites = I(:);
%% 1. Data term only
labels = initMRF(sites);
MAP_labels1 = do_greedyMAP(labels);

%% 2. with Pairwise term and ICM
MAP_labels2 = icm(labels);

%% 3. with Simulated Annealing
MAP_labels3 = sim_annealing(labels)

%% 4. CRF
 
