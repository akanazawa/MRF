%%%%%%%%%%%%%%%%%%%%
% MRF Problem set 3 CMSC 828 Spring '12
% 
%
% Angjoo Kanazawa 4/30/'12
%%%%%%%%%%%%%%%%%%%%
addpath('utils');
I = double(imread('TestImage.jpg'));
sites = initSites(I);

%% 1. Data term only
[U, initLabels] = mapUnary(I(:));

%% plot
img1 = reshape(initLabels, size(I));
sfigure; imshow(img1, []);

%% 2. with Pairwise term and ICM
K = .897; 
beta = .95;

icmLabels = icm(sites, initLabels, I, K, beta, []);
img2 = reshape(icmLabels, size(I));
sfigure; imshow(img2, []);


%% 3. with Simulated Annealing
saLabels = sim_annealing(sites, initLabels, K, beta);
img3 = reshape(saLabels, size(I));
sfigure; imshow(img3, []);

%% 4. CRF
sigma = 1.5;
I2 = gaussfilt(I, sigma);
[Ix Iy] = derivative5(I2, 'x', 'y');        
CRF.Ix = Ix; CRF.Iy = Iy;
CRF.Mx = max(max(abs(Ix)));
CRF.My = max(max(abs(Iy)));

crfLabels = icm(sites, initLabels, I, K, beta, CRF);
img4 = reshape(saLabels, size(I));
sfigure; imshow(img4, []);
